#!/usr/bin/env python3
from pathlib import Path
import pandas as pd

ROOT = Path(__import__("os").environ.get("HFD_PROJECT_ROOT", ".")).resolve()

P92 = ROOT / "00_metadata/human_validation/analysis_outputs/gene_meta_v1/human_MASLD_gene_meta_summary_v1.txt"
P93 = ROOT / "00_metadata/human_validation/analysis_outputs/mouse_human_concordance_v1_6/mouse_human_reproducibility_concordance_summary_v1_6.txt"
P94 = ROOT / "00_metadata/human_validation/analysis_outputs/mouse_human_concordance_followup_v1_1/mouse_human_concordance_shape_summary_v1_1.txt"
P95 = ROOT / "00_metadata/human_validation/analysis_outputs/pathway_concordance_v1_1/mouse_human_Hallmark_concordance_summary_v1_1.txt"

PATHJOIN = ROOT / "07_human_validation/pathway_concordance_v1_1/mouse_HFD_human_MASLD_Hallmark_concordance_v1_1.tsv"
GENEJOIN = ROOT / "07_human_validation/mouse_human_concordance_v1_6/mouse_HFD_human_MASLD_one2one_ortholog_concordance_v1_6.tsv.gz"

OUTDIR = ROOT / "07_human_validation/final_summary_v1"
METADIR = ROOT / "00_metadata/human_validation/analysis_outputs/final_summary_v1"
OUTDIR.mkdir(parents=True, exist_ok=True)
METADIR.mkdir(parents=True, exist_ok=True)

OUT_PATH = OUTDIR / "mouse_human_Hallmark_joint_significance_classes_v1.tsv"
OUT_Q4 = OUTDIR / "mouse_high_reproducibility_human_concordant_genes_top_v1.tsv.gz"
OUT_SUM = METADIR / "human_validation_integrated_summary_v1.txt"

for x in [P92,P93,P94,P95,PATHJOIN,GENEJOIN]:
    if not x.exists():
        raise FileNotFoundError(f"Missing required input: {x}")

def parse_summary(path):
    d = {}
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("=") or "=" not in line:
            continue
        k,v = line.split("=",1)
        d[k.strip()] = v.strip()
    return d

s92 = parse_summary(P92)
s93 = parse_summary(P93)
s94 = parse_summary(P94)
s95 = parse_summary(P95)

pw = pd.read_csv(PATHJOIN, sep="\t")

required_pw = {
    "pathway","mouse_NES","mouse_FDR","human_RE_NES","human_RE_FDR",
    "same_direction_RE","both_FDR05_RE","both_FDR05_same_direction_RE"
}
missing = required_pw - set(pw.columns)
if missing:
    raise RuntimeError(f"Pathway join missing columns: {sorted(missing)}")

pw["joint_class"] = "NEITHER_FDR05"
pw.loc[(pw["mouse_FDR"] < 0.05) & (pw["human_RE_FDR"] >= 0.05), "joint_class"] = "MOUSE_ONLY_FDR05"
pw.loc[(pw["mouse_FDR"] >= 0.05) & (pw["human_RE_FDR"] < 0.05), "joint_class"] = "HUMAN_ONLY_FDR05"
pw.loc[(pw["mouse_FDR"] < 0.05) & (pw["human_RE_FDR"] < 0.05) & (pw["same_direction_RE"] == True), "joint_class"] = "BOTH_FDR05_SAME_DIRECTION"
pw.loc[(pw["mouse_FDR"] < 0.05) & (pw["human_RE_FDR"] < 0.05) & (pw["same_direction_RE"] == False), "joint_class"] = "BOTH_FDR05_OPPOSITE_DIRECTION"

pw["abs_human_RE_NES"] = pw["human_RE_NES"].abs()
pw = pw.sort_values(["joint_class","human_RE_FDR","abs_human_RE_NES"], ascending=[True,True,False])
pw = pw.drop(columns=["abs_human_RE_NES"])
pw.to_csv(OUT_PATH, sep="\t", index=False)

# Gene-level top-Q4, human-concordant descriptive table.
g = pd.read_csv(GENEJOIN, sep="\t", compression="gzip")
required_g = {
    "mouse_ensembl","mouse_symbol","human_ENTREZID","human_symbol",
    "mouse_log2FC","mouse_direction_consistency","mouse_nGSE",
    "human_RE_log2FC","human_RE_FDR","concordant_primary",
    "human_aligned_effect"
}
missing = required_g - set(g.columns)
if missing:
    raise RuntimeError(f"Gene join missing columns: {sorted(missing)}")

q4_cut = g["mouse_direction_consistency"].quantile(0.75)
q4 = g[(g["mouse_direction_consistency"] >= q4_cut) & (g["concordant_primary"] == True)].copy()
q4["human_FDR05"] = q4["human_RE_FDR"] < 0.05
q4 = q4.sort_values(
    ["human_FDR05","mouse_direction_consistency","human_RE_FDR","human_aligned_effect"],
    ascending=[False,False,True,False]
)
q4.to_csv(OUT_Q4, sep="\t", index=False, compression="gzip")

joint_counts = pw["joint_class"].value_counts().to_dict()

both_same = pw[pw["joint_class"]=="BOTH_FDR05_SAME_DIRECTION"]["pathway"].tolist()
both_opp = pw[pw["joint_class"]=="BOTH_FDR05_OPPOSITE_DIRECTION"]["pathway"].tolist()

summary = f"""===== INTEGRATED HUMAN VALIDATION SUMMARY v1 =====

STATUS=PASS
ROLE=FINAL_TRANSLATIONAL_RESULT_FREEZE
UPSTREAM_HUMAN_META=92_v1
UPSTREAM_GENE_CONCORDANCE=93_v1_6
UPSTREAM_SHAPE_CHARACTERIZATION=94_v1_1
UPSTREAM_PATHWAY_CONCORDANCE=95_v1_1
PRIMARY_TRANSLATIONAL_INFERENCE_CHANGED=NO

HUMAN_COHORT_META:
Primary_common_genes={s92.get('Primary_common_genes','NA')}
Primary_RE_FDR_lt_0.05={s92.get('Primary_RE_FDR_lt_0.05','NA')}
Primary_median_I2_percent={s92.get('Primary_median_I2_percent','NA')}
Primary_vs_strict4_RE_effect_Pearson_r={s92.get('Primary_vs_strict4_RE_effect_Pearson_r','NA')}
Primary_vs_strict4_same_direction_fraction={s92.get('Primary_vs_strict4_same_direction_fraction','NA')}

GENE_LEVEL_MOUSE_HUMAN:
Final_joined_translational_genes={s93.get('Final_joined_translational_genes','NA')}
Overall_human_same_direction_fraction={s93.get('Overall_human_same_direction_fraction','NA')}
Primary_OR_per_10pp={s93.get('Primary_OR_per_10pp','NA')}
Primary_95CI={s93.get('Primary_95CI','NA')}
Primary_P={s93.get('Primary_P','NA')}
Adjusted_OR_per_10pp={s93.get('Adjusted_OR_per_10pp','NA')}
Adjusted_P={s93.get('Adjusted_P','NA')}
Coverage_ge20_OR_per_10pp={s93.get('Coverage_ge20_OR_per_10pp','NA')}
Coverage_ge20_P={s93.get('Coverage_ge20_P','NA')}
Human_FE_sensitivity_OR_per_10pp={s93.get('Human_FE_sensitivity_OR_per_10pp','NA')}
Strict4_sensitivity_OR_per_10pp={s93.get('Strict4_sensitivity_OR_per_10pp','NA')}
Mouse_vs_human_effect_Pearson_r={s93.get('Mouse_vs_human_effect_Pearson_r','NA')}
Mouse_vs_human_effect_Spearman_rho={s93.get('Mouse_vs_human_effect_Spearman_rho','NA')}

SECONDARY_SHAPE:
Q4_vs_Q1_unadjusted_OR={s94.get('Q4_vs_Q1_unadjusted_OR','NA')}
Q4_vs_Q1_unadjusted_95CI={s94.get('Q4_vs_Q1_unadjusted_95CI','NA')}
Q4_vs_Q1_unadjusted_P={s94.get('Q4_vs_Q1_unadjusted_P','NA')}
Q4_vs_lower75_adjusted_OR={s94.get('Q4_vs_lower75_adjusted_OR','NA')}
Q4_vs_lower75_adjusted_P={s94.get('Q4_vs_lower75_adjusted_P','NA')}
PI_excludes_zero_N={s94.get('PI_excludes_zero_N','NA')}
PI_excludes_zero_adjusted_OR={s94.get('PI_excludes_zero_adjusted_OR','NA')}
PI_excludes_zero_adjusted_P={s94.get('PI_excludes_zero_adjusted_P','NA')}

PATHWAY_LEVEL_MOUSE_HUMAN:
Hallmark_pathways_joined={s95.get('Hallmark_pathways_joined','NA')}
Human_RE_Hallmark_FDR_lt_0.05={s95.get('Human_RE_Hallmark_FDR_lt_0.05','NA')}
Mouse_vs_human_RE_NES_Pearson_r={s95.get('Mouse_vs_human_RE_NES_Pearson_r','NA')}
Mouse_vs_human_RE_NES_Pearson_95CI={s95.get('Mouse_vs_human_RE_NES_Pearson_95CI','NA')}
Mouse_vs_human_RE_NES_Spearman_rho={s95.get('Mouse_vs_human_RE_NES_Spearman_rho','NA')}
Mouse_vs_human_RE_same_direction_fraction={s95.get('Mouse_vs_human_RE_same_direction_fraction','NA')}
Mouse_and_human_RE_both_FDR05={s95.get('Mouse_and_human_RE_both_FDR05','NA')}
Mouse_and_human_RE_both_FDR05_same_direction={s95.get('Mouse_and_human_RE_both_FDR05_same_direction','NA')}

PATHWAY_JOINT_CLASSES:
BOTH_FDR05_SAME_DIRECTION={joint_counts.get('BOTH_FDR05_SAME_DIRECTION',0)}
BOTH_FDR05_OPPOSITE_DIRECTION={joint_counts.get('BOTH_FDR05_OPPOSITE_DIRECTION',0)}
MOUSE_ONLY_FDR05={joint_counts.get('MOUSE_ONLY_FDR05',0)}
HUMAN_ONLY_FDR05={joint_counts.get('HUMAN_ONLY_FDR05',0)}
NEITHER_FDR05={joint_counts.get('NEITHER_FDR05',0)}

BOTH_FDR05_SAME_DIRECTION_PATHWAYS={';'.join(both_same)}
BOTH_FDR05_OPPOSITE_DIRECTION_PATHWAYS={';'.join(both_opp)}

INTERPRETATION:
1=Overall_gene_level_mouse_human_effect_concordance_is_weak
2=Higher_mouse_cross_experiment_directional_reproducibility_is_associated_with_modestly_higher_human_directional_concordance
3=Association_is_enriched_in_the_highest_reproducibility_quartile_but_primary_continuous_model_remains_prespecified_inference
4=Mouse_prediction_interval_exclusion_does_not_show_additional_translational_predictive_value_in_this_sparse_subset
5=Hallmark_level_global_concordance_is_not_supported
6=Translational_signal_is_selective_not_global
7=CAUSAL_INTERPRETATION_NO

OUTPUT_PATHWAY_CLASSES={OUT_PATH.relative_to(ROOT)}
OUTPUT_Q4_CONCORDANT_GENES={OUT_Q4.relative_to(ROOT)}
"""

OUT_SUM.write_text(summary)
print(summary)

