"""build_data_source_table.py

Builds the data-source table (STROBE-MR items 4 and 10) from data/raw/MANIFEST.tsv. No values are typed by hand:
sample sizes, ancestry, accession and UKB inclusion are copied from the manifest; only the trait label and the role in
this study are assigned here.
Usage (repository root): python3 manuscript/build_data_source_table.py
Output: manuscript/tables/table_data_sources.tsv and .md
"""
import pandas as pd

M = pd.read_csv("data/raw/MANIFEST.tsv", sep="\t", dtype=str)
USED = [  # (manifest dataset, trait, role in this study)
    ("ea_okbay2022", "Educational attainment (EA4)", "Exposure, primary"),
    ("income_kweon2025", "Income factor", "Exposure, primary"),
    ("tdi_panukb", "Townsend deprivation index", "Exposure, primary (all ancestries; EUR-only file as sensitivity)"),
    ("occupation_akimova2024", "Occupational status (CAMSIS)", "Exposure, primary"),
    ("ea_lee2018", "Educational attainment (EA3)", "Exposure, sensitivity"),
    ("income_hill2019", "Household income", "Exposure, sensitivity; Li D-comparable model"),
    ("occupation_akimova2024_isei", "Occupational status (ISEI)", "Exposure, Le-comparable model only"),
    ("ses_gses_xia2025", "General SES factor (gSES)", "Exposure, secondary"),
    ("cad_aragam2022_eur", "Coronary artery disease", "Outcome, primary"),
    ("cad_aragam2022_multi", "Coronary artery disease", "Outcome, comparator"),
    ("cad_mvp_phe411_eur", "Ischaemic heart disease (phecode 411)", "Outcome, comparator without UKB"),
    ("pad_mvp_phe440_2_eur", "Peripheral artery disease (phecode 440.2)", "Outcome, primary"),
    ("pad_finngen_r12_athscle", "Atherosclerosis (I9_ATHSCLE)", "Outcome, comparator"),
    ("pad_vanzuydam2021", "Peripheral artery disease", "Outcome, exploratory (rounded data)"),
    ("pad_mvp_phe440_2_afr", "Peripheral artery disease (phecode 440.2)", "Outcome, exploratory transportability"),
    ("hf_hermes2_2025", "Heart failure", "Outcome, primary"),
    ("hf_mvp_phe428_eur", "Heart failure (phecode 428)", "Outcome, comparator without UKB"),
    ("stroke_gigastroke_any_multi", "Any stroke", "Outcome, primary"),
    ("stroke_gigastroke_any_eur", "Any stroke", "Outcome, comparator"),
    ("stroke_mvp_phe433_21_eur", "Cerebral artery occlusion with infarction (phecode 433.21)", "Outcome, comparator without UKB"),
    ("bmi_giant_yengo2018", "Body mass index", "Mediator (exploratory)"),
    ("smoking_wootton2020", "Lifetime smoking index", "Mediator (exploratory)"),
    ("ldl_glgc_graham2021", "LDL cholesterol", "Mediator (exploratory)"),
    ("hdl_glgc_graham2021", "HDL cholesterol", "Mediator (exploratory)"),
    ("tg_glgc_graham2021", "Triglycerides", "Mediator (exploratory)"),
    ("sbp_keaton2024", "Systolic blood pressure", "Mediator (exploratory)"),
    ("t2d_diamante_mahajan2022", "Type 2 diabetes", "Mediator (exploratory)"),
]
rows = []
for ds, trait, role in USED:
    r = M[M.dataset == ds]
    assert len(r) == 1, ds
    r = r.iloc[0]
    rows.append({"trait": trait, "role": role, "source": r.citation, "accession": r.accession, "ancestry": r.ancestry,
                 "n_total": r.total_n, "cases": r.cases, "controls": r.controls, "includes_ukb": r.includes_ukb,
                 "manifest_dataset": ds})
T = pd.DataFrame(rows).fillna("")
# HDL and TG come from the same GLGC 2021 trans-ancestry meta-analysis as LDL; use the manifest's LDL ancestry breakdown.
ldl_anc = T.loc[T.manifest_dataset == "ldl_glgc_graham2021", "ancestry"].iloc[0]
T.loc[T.manifest_dataset.isin(["hdl_glgc_graham2021", "tg_glgc_graham2021"]), "ancestry"] = ldl_anc
T.to_csv("manuscript/tables/table_data_sources.tsv", sep="\t", index=False)
cols = ["trait", "role", "source", "ancestry", "n_total", "cases", "controls", "includes_ukb"]
with open("manuscript/tables/table_data_sources.md", "w") as f:
    f.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
    for _, r in T.iterrows():
        f.write("| " + " | ".join(str(r[c]).replace("|", "/") for c in cols) + " |\n")
print(len(T))

# GWAS covariate table (STROBE-MR items 4 and 6): source-checked entries compiled 2026-10-08 in
# manuscript/verification/covariate_sources/ (one row per manifest dataset; 'undocumented in ...' marks gaps).
C = pd.concat([pd.read_csv(f"manuscript/verification/covariate_sources/{f}", sep="\t", dtype=str) for f in
               ["gwas_covariates_exposures_mediators.tsv", "gwas_covariates_outcomes.tsv"]])
assert set(C.manifest_dataset) == set(T.manifest_dataset), set(C.manifest_dataset) ^ set(T.manifest_dataset)
C = T[["manifest_dataset", "trait", "role"]].merge(C.drop(columns=[c for c in ["trait"] if c in C.columns]), on="manifest_dataset")
C.to_csv("manuscript/tables/table_gwas_covariates.tsv", sep="\t", index=False)
print(len(C), C.status.value_counts().to_dict())
