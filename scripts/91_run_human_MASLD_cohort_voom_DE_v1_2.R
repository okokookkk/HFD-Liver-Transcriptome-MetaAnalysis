#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  req <- c("data.table","edgeR","limma","AnnotationDbi","org.Hs.eg.db")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly=TRUE)]
  if (length(miss)) stop("Missing required R packages: ", paste(miss, collapse=", "))
})
library(data.table); library(edgeR); library(limma); library(AnnotationDbi); library(org.Hs.eg.db)

ROOT <- normalizePath(Sys.getenv("HFD_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
GROUPS <- file.path(ROOT,"00_metadata/human_validation/frozen_groups_v1_2/human_MASLD_primary_sample_groups_frozen_v1_2.tsv")
G130S <- file.path(ROOT,"00_metadata/human_validation/frozen_groups_v1_2/GSE130970_strict4_control_sensitivity_groups_v1_2.tsv")
MAPF <- file.path(ROOT,"00_metadata/human_validation/frozen_groups_v1_2/human_MASLD_count_sample_mapping_audit_v1_2.tsv")
G126 <- file.path(ROOT,"00_metadata/human_validation/supplementary_v1/GSE126848/GSE126848_Gene_counts_raw.txt.gz")
G130 <- file.path(ROOT,"00_metadata/human_validation/supplementary_v1/GSE130970/GSE130970_all_sample_salmon_tximport_counts_entrez_gene_ID.csv.gz")
G135T <- file.path(ROOT,"00_metadata/human_validation/supplementary_v1/GSE135251/GSE135251_RAW.tar")

OUTDIR <- file.path(ROOT,"07_human_validation/cohort_DE_v1_2")
METADIR <- file.path(ROOT,"00_metadata/human_validation/analysis_outputs/cohort_DE_v1_2")
dir.create(OUTDIR,recursive=TRUE,showWarnings=FALSE); dir.create(METADIR,recursive=TRUE,showWarnings=FALSE)

for (q in c(GROUPS,G130S,MAPF,G126,G130,G135T)) if(!file.exists(q)) stop("Missing input: ",q)

groups <- fread(GROUPS); g130s <- fread(G130S); mp <- fread(MAPF)

map_ens <- function(ids){
  e <- sub("\\.[0-9]+$","",as.character(ids))
  z <- as.data.table(AnnotationDbi::select(org.Hs.eg.db, keys=unique(e), keytype="ENSEMBL", columns="ENTREZID"))
  z <- unique(z[!is.na(ENTREZID) & nzchar(ENTREZID),.(ENSEMBL,ENTREZID)])
  c1 <- z[,.(n=uniqueN(ENTREZID)),by=ENSEMBL]
  z <- z[ENSEMBL %in% c1[n==1,ENSEMBL]]
  list(e=e,map=z,n_input=uniqueN(e),n_map=uniqueN(z$ENSEMBL),n_amb=sum(c1$n>1))
}
collapse_ens <- function(mat,ids){
  mm <- map_ens(ids); ix <- match(mm$e,mm$map$ENSEMBL); keep <- !is.na(ix)
  entrez <- mm$map$ENTREZID[ix[keep]]
  m <- mat[keep,,drop=FALSE]; storage.mode(m)<-"double"
  list(counts=rowsum(m,entrez,reorder=FALSE),
       audit=c(source=nrow(mat),unique_ensembl=mm$n_input,mapped_ensembl=mm$n_map,
               ambiguous_ensembl=mm$n_amb,entrez=nrow(rowsum(m,entrez,reorder=FALSE))))
}
collapse_entrez <- function(mat,ids){
  ids <- trimws(as.character(ids)); keep <- grepl("^[0-9]+$",ids)
  m <- mat[keep,,drop=FALSE]; storage.mode(m)<-"double"
  cmat <- rowsum(m,ids[keep],reorder=FALSE)
  list(counts=cmat,audit=c(source=nrow(mat),numeric_source=sum(keep),dropped=sum(!keep),entrez=nrow(cmat)))
}

canonical_count_id <- function(x) {
  x <- trimws(as.character(x))
  out <- x
  isnum <- grepl("^[0-9]+$", x)
  if (any(isnum)) {
    # numeric equivalence only for pure digit IDs; preserve all other identifiers.
    z <- sub("^0+", "", x[isnum])
    z[z == ""] <- "0"
    out[isnum] <- z
  }
  out
}

map_count_headers_to_gsm <- function(raw_headers, z, gse, allow_numeric=FALSE) {
  raw_headers <- as.character(raw_headers)
  frozen_ids <- as.character(z$Count_sample_identifier)
  frozen_gsm <- as.character(z$Mapped_GSM)

  if (any(is.na(frozen_gsm)) || any(frozen_gsm == ""))
    stop(gse, ": frozen mapping contains empty GSM")

  # First try exact string identity.
  ix <- match(raw_headers, frozen_ids)
  exact_n <- sum(!is.na(ix))

  if (allow_numeric && anyNA(ix)) {
    raw_key <- canonical_count_id(raw_headers)
    frozen_key <- canonical_count_id(frozen_ids)

    if (anyDuplicated(frozen_key))
      stop(gse, ": canonical frozen count IDs are not unique")

    miss <- which(is.na(ix))
    ix2 <- match(raw_key[miss], frozen_key)
    ix[miss] <- ix2
  }

  if (anyNA(ix)) {
    bad <- raw_headers[is.na(ix)]
    stop(
      gse, ": frozen mapping mismatch after permitted ID normalization; ",
      "unmapped=", length(bad), "; examples=", paste(head(bad, 10L), collapse=",")
    )
  }

  mapped_gsm <- frozen_gsm[ix]

  if (anyDuplicated(mapped_gsm))
    stop(gse, ": duplicate GSM after frozen count-ID mapping")

  list(
    gsm = mapped_gsm,
    n_exact = exact_n,
    n_normalized = length(raw_headers) - exact_n
  )
}

load126 <- function(){
  x<-fread(G126, header=TRUE, check.names=FALSE)
  if (ncol(x) != 58L) stop("GSE126848 expected 58 columns after forced header parse, got ", ncol(x))
  ids<-x[[1]]
  raw_headers <- names(x)[-1]
  if (all(grepl("^V[0-9]+$", raw_headers)))
    stop("GSE126848 header parse still failed: generic V2... names remain")
  m<-as.matrix(x[,-1,with=FALSE])
  colnames(m)<-raw_headers
  z<-mp[GSE=="GSE126848"]
  mm <- map_count_headers_to_gsm(
    colnames(m), z, "GSE126848", allow_numeric=TRUE
  )
  colnames(m) <- mm$gsm

  frozen_gsm <- groups[GSE=="GSE126848", GSM]
  if (!setequal(colnames(m), frozen_gsm))
    stop("GSE126848 mapped count GSM set != frozen biological GSM set")

  cc <- collapse_ens(m,ids)
  cc$mapping_audit <- c(
    count_columns=ncol(m),
    exact_header_matches=mm$n_exact,
    numeric_equivalence_matches=mm$n_normalized,
    mapped_unique_GSM=uniqueN(colnames(m))
  )
  cc
}
load130 <- function(){
  x<-fread(G130, header=TRUE, check.names=FALSE)
  if (ncol(x) != 79L) stop("GSE130970 expected 79 columns after forced header parse, got ", ncol(x))
  ids<-x[[1]]
  raw_headers <- names(x)[-1]
  m<-as.matrix(x[,-1,with=FALSE])
  colnames(m)<-raw_headers
  z<-mp[GSE=="GSE130970"]
  mm <- map_count_headers_to_gsm(
    colnames(m), z, "GSE130970", allow_numeric=FALSE
  )
  colnames(m) <- mm$gsm

  frozen_gsm <- groups[GSE=="GSE130970", GSM]
  if (!setequal(colnames(m), frozen_gsm))
    stop("GSE130970 mapped count GSM set != frozen biological GSM set")

  cc <- collapse_entrez(m,ids)
  cc$mapping_audit <- c(
    count_columns=ncol(m),
    exact_header_matches=mm$n_exact,
    normalized_matches=mm$n_normalized,
    mapped_unique_GSM=uniqueN(colnames(m))
  )
  cc
}
load135 <- function(){
  td<-tempfile("g135_"); dir.create(td); on.exit(unlink(td,recursive=TRUE,force=TRUE),add=TRUE)
  untar(G135T,exdir=td)
  allfs <- list.files(td,recursive=TRUE,full.names=TRUE)
  zmap <- mp[GSE=="GSE135251"]
  frozen_names <- as.character(zmap$Count_sample_identifier)
  fs <- allfs[basename(allfs) %in% frozen_names]
  if(length(fs)!=216)
    stop("Expected 216 frozen GSE135251 files, got ",length(fs))

  # Reorder files to frozen mapping order, then take GSM directly from the frozen map.
  j <- match(frozen_names, basename(fs))
  if(anyNA(j)) stop("GSE135251 frozen filename mapping mismatch")
  fs <- fs[j]
  gsm <- as.character(zmap$Mapped_GSM)

  if(any(is.na(gsm)) || any(gsm=="") || anyDuplicated(gsm))
    stop("GSE135251 invalid frozen GSM mapping")
  z0<-fread(fs[1],header=FALSE); ids0<-as.character(z0[[1]]); v0<-suppressWarnings(as.numeric(z0[[2]]))
  keep0<-!grepl("^__",ids0) & is.finite(v0); ids<-ids0[keep0]
  m<-matrix(NA_real_,nrow=length(ids),ncol=length(fs),dimnames=list(ids,gsm))
  for(i in seq_along(fs)){
    z<-fread(fs[i],header=FALSE); f<-as.character(z[[1]]); v<-suppressWarnings(as.numeric(z[[2]]))
    keep<-!grepl("^__",f) & is.finite(v); f<-f[keep]; v<-v[keep]
    if(!identical(f,ids)){ j<-match(ids,f); if(anyNA(j)) stop("GSE135251 feature mismatch"); v<-v[j] }
    m[,i]<-v
  }
  cc <- collapse_ens(m,ids)
  cc$mapping_audit <- c(
    count_files=ncol(m),
    frozen_filename_matches=length(fs),
    mapped_unique_GSM=uniqueN(colnames(m))
  )
  cc
}

align <- function(cmat,gse,gt){
  g<-gt[GSE==gse,.(GSM,Human_group)]
  if(!setequal(colnames(cmat),g$GSM)) stop(gse,": count/group GSM mismatch")
  g<-g[match(colnames(cmat),GSM)]
  list(counts=cmat,group=factor(g$Human_group,levels=c("CONTROL","CASE")))
}
run_de <- function(cmat,grp,label){
  y<-DGEList(cmat,group=grp)
  design<-model.matrix(~0+grp); colnames(design)<-c("CONTROL","CASE")
  keep<-filterByExpr(y,design=design); y<-y[keep,,keep.lib.sizes=FALSE]
  y<-calcNormFactors(y,method="TMM"); v<-voom(y,design,plot=FALSE)
  fit<-lmFit(v,design); ct<-makeContrasts(CASE_minus_CONTROL=CASE-CONTROL,levels=design)
  fit<-eBayes(contrasts.fit(fit,ct))
  b<-as.numeric(fit$coefficients[,1]); se<-as.numeric(fit$stdev.unscaled[,1]*sqrt(fit$s2.post))
  pv<-as.numeric(fit$p.value[,1])
  out<-data.table(ENTREZID=rownames(fit$coefficients),human_log2FC=b,human_SE=se,
                  human_sampling_variance=se^2,human_t=as.numeric(fit$t[,1]),
                  human_p_value=pv,human_FDR_BH=p.adjust(pv,"BH"),
                  AveExpr_log2CPM=as.numeric(fit$Amean),N_case=sum(grp=="CASE"),
                  N_control=sum(grp=="CONTROL"),Analysis=label)
  list(out=out,n_test=nrow(out),n_sig=sum(out$human_FDR_BH<.05),
       n_up=sum(out$human_FDR_BH<.05 & out$human_log2FC>0),
       n_dn=sum(out$human_FDR_BH<.05 & out$human_log2FC<0))
}

cat("===== HUMAN MASLD COHORT VOOM DE v1.2 =====\n")
a126<-load126(); x126<-align(a126$counts,"GSE126848",groups); d126<-run_de(x126$counts,x126$group,"GSE126848_PRIMARY_31v26")
a130<-load130(); x130<-align(a130$counts,"GSE130970",groups); d130<-run_de(x130$counts,x130$group,"GSE130970_PRIMARY_72v6")
x130s<-align(a130$counts,"GSE130970",g130s); d130s<-run_de(x130s$counts,x130s$group,"GSE130970_STRICT4_74v4")
a135<-load135(); x135<-align(a135$counts,"GSE135251",groups); d135<-run_de(x135$counts,x135$group,"GSE135251_PRIMARY_206v10")

outs<-list(
 GSE126848=file.path(OUTDIR,"GSE126848_MASLD_vs_control_voom_DE_v1_2.tsv.gz"),
 GSE130970=file.path(OUTDIR,"GSE130970_MASLD_vs_control6_voom_DE_v1_2.tsv.gz"),
 GSE130970_strict4=file.path(OUTDIR,"GSE130970_MASLD_vs_strict_control4_voom_DE_v1_2.tsv.gz"),
 GSE135251=file.path(OUTDIR,"GSE135251_MASLD_vs_control_voom_DE_v1_2.tsv.gz")
)
fwrite(d126$out,outs$GSE126848,sep="\t",compress="gzip")
fwrite(d130$out,outs$GSE130970,sep="\t",compress="gzip")
fwrite(d130s$out,outs$GSE130970_strict4,sep="\t",compress="gzip")
fwrite(d135$out,outs$GSE135251,sep="\t",compress="gzip")

fmt <- function(x) paste(names(x),x,sep="=",collapse=";")
sumfile<-file.path(METADIR,"human_MASLD_cohort_voom_DE_summary_v1_2.txt")
summary<-c(
"===== HUMAN MASLD COHORT VOOM DE SUMMARY v1.2 =====","",
"STATUS=PASS",
"SUPERSEDES=v1_and_v1_1_loader_header_detection",
"GSE126848_FREAD_HEADER_FORCED=TRUE",
"GSE126848_FREAD_CHECK_NAMES=FALSE",
"FROZEN_HUMAN_GROUPS=00_metadata/human_validation/frozen_groups_v1_2",
paste0("edgeR_version=",packageVersion("edgeR")),
paste0("limma_version=",packageVersion("limma")),
"PRIMARY_DESIGN=group_only_CASE_minus_CONTROL",
"NORMALIZATION=TMM","MODEL=limma_voom_empirical_Bayes","FILTER=edgeR_filterByExpr",
"GSE130970_COUNTS_ROUNDED=NO","CAUSAL_INTERPRETATION=NO","",
paste0("GSE126848_mapping_audit=",fmt(a126$mapping_audit)),
paste0("GSE126848_feature_audit=",fmt(a126$audit)),
paste0("GSE126848_genes_tested=",d126$n_test),
paste0("GSE126848_FDR_lt_0.05=",d126$n_sig),
paste0("GSE126848_FDR05_up=",d126$n_up),
paste0("GSE126848_FDR05_down=",d126$n_dn),"",
paste0("GSE130970_mapping_audit=",fmt(a130$mapping_audit)),
paste0("GSE130970_feature_audit=",fmt(a130$audit)),
paste0("GSE130970_primary_genes_tested=",d130$n_test),
paste0("GSE130970_primary_FDR_lt_0.05=",d130$n_sig),
paste0("GSE130970_strict4_genes_tested=",d130s$n_test),
paste0("GSE130970_strict4_FDR_lt_0.05=",d130s$n_sig),"",
paste0("GSE135251_mapping_audit=",fmt(a135$mapping_audit)),
paste0("GSE135251_feature_audit=",fmt(a135$audit)),
paste0("GSE135251_genes_tested=",d135$n_test),
paste0("GSE135251_FDR_lt_0.05=",d135$n_sig),
paste0("GSE135251_FDR05_up=",d135$n_up),
paste0("GSE135251_FDR05_down=",d135$n_dn),"",
"Human_meta_analysis_performed=NO","Mouse_human_concordance_performed=NO"
)
writeLines(summary,sumfile); cat(paste(summary,collapse="\n"),"\n")

