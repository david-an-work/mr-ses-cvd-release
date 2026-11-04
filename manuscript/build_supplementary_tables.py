"""build_supplementary_tables.py

Assembles the supplementary tables workbook from the saved results tables (no computation of new estimates).
Usage (from the repository root): python3 manuscript/build_supplementary_tables.py
Output: manuscript/supplementary_tables.xlsx (one sheet per table, plus a contents sheet)
"""
import os
import pandas as pd

T = "results/tables"
OUT = "manuscript/supplementary_tables.xlsx"

def rd(name, **kw):
    return pd.read_csv(os.path.join(T, name), sep="\t", **kw)

sheets = []  # (sheet name, title, source, dataframe); sheets are numbered S1, S2, ... in list order at write time
sheets.append(("data_sources", "GWAS data sources: role, ancestry, sample size and UK Biobank inclusion (from data/raw/MANIFEST.tsv)", "manuscript/tables/table_data_sources.tsv", pd.read_csv("manuscript/tables/table_data_sources.tsv", sep="\t", dtype=str)))

sheets.append(("gwas_covariates", "GWAS ancestry handling, phenotype preprocessing and association-model covariates per source (source-checked; gaps marked undocumented)", "manuscript/tables/table_gwas_covariates.tsv", pd.read_csv("manuscript/tables/table_gwas_covariates.tsv", sep="\t", dtype=str)))
inst = pd.read_csv("data/instruments/INSTRUMENT_SUMMARY.tsv", sep="\t")
sheets.append(("instruments", "Instrument selection per exposure and mediator", "data/instruments/INSTRUMENT_SUMMARY.tsv", inst))
sheets.append(("genetic_corr", "Cross-trait LDSC genetic correlations among the SES measures (gSES correlations in Supplementary Table S20)", "feasibility_rg_ldsc.tsv", rd("feasibility_rg_ldsc.tsv")))
sheets.append(("harmonization", "Instrument retention per harmonized cell (primary analyses, no proxies). Includes harmonized sets from the registered mediator-strength feasibility check (SES measures other than EA4 with mediators, and full mediator panels) that are not reported in the main text; the manuscript uses the univariable cells, the 11 SES-versus-SES sets and the 7 EA4 + mediator sets", "harmonization_retention.tsv", rd("harmonization_retention.tsv")))

est = rd("mr_univariable_estimates.tsv")
sheets.append(("univariable_all", "Univariable MR, all methods, all 96 cells (OR per SD; gSES per released factor unit; Hill income per released unit, not SD)", "mr_univariable_estimates.tsv", est))
grid = rd("mr_primary_grid.tsv")
sheets.append(("primary_grid", "Primary 16-cell grid (IVW, multiplicative random effects) with Bonferroni-level CIs", "mr_primary_grid.tsv", grid))
dg = rd("mr_univariable_diagnostics.tsv")
sheets.append(("diagnostics", "Heterogeneity, MR-Egger intercept, MR-PRESSO, leave-one-out and MR-Steiger summaries", "mr_univariable_diagnostics.tsv", dg))
sheets.append(("sensitivity_variants", "IVW sensitivity variants (proxies, MHC excluded, LD-pruned, MVP r2 >= 0.8, Steiger-filtered, van Zuydam SE)", "mr_univariable_variants.tsv", rd("mr_univariable_variants.tsv")))
sheets.append(("appraisal", "SES-versus-SES MVMR: conditional F at four covariance scales and IVW-MVMR estimates (all weak-instrument-limited)", "mvmr_appraisal.tsv", rd("mvmr_appraisal.tsv")))
med = rd("mediation_education.tsv")
sheets.append(("mediation", "Education single-mediator models: indirect, difference and ratio summaries (exploratory; unadjusted intervals)", "mediation_education.tsv", med))
sheets.append(("mediation_condF", "Conditional F for the single-mediator models (strength gate)", "mediation_education_condF.tsv", rd("mediation_education_condF.tsv")))
sheets.append(("overlap_scenarios", "Illustrative Burgess-type overlap scenarios (phi_UKB x 1/F; not bias estimates)", "overlap_scenarios.tsv", rd("overlap_scenarios.tsv")))
if os.path.exists(os.path.join(T, "overlap_ldsc_exposure_outcome.tsv")):
    sheets.append(("overlap_ldsc", "Exposure x outcome cross-trait LDSC intercepts (descriptive)", "overlap_ldsc_exposure_outcome.tsv", rd("overlap_ldsc_exposure_outcome.tsv")))
pv = os.path.join(T, "mr_presso_outlier_verification.tsv")
if os.path.exists(pv):
    sheets.append(("presso_verification", "Re-run of MR-PRESSO for cells with corrected outlier counts", "mr_presso_outlier_verification.tsv", rd("mr_presso_outlier_verification.tsv")))

sheets.append(("power_univariable", "Pre-estimate planning benchmarks: expected IVW SE and OR per SD detectable with 80% power at alpha 0.003125 for every cell (Hill income and gSES: per released unit; computed from outcome SEs and discovery exposure effects only; not post hoc power)", "outcome_precision_univariable.tsv", rd("outcome_precision_univariable.tsv")))
sheets.append(("precision_mvmr", "Pre-estimate expected fixed-effect IVW-MVMR SEs (planning values) for the 8 models for which planning SEs were computed before estimation (labels as in the plan: 4-SES = A1, 4-SES (Hill) = A5, EA4+Kweon = A2, EA4+TDI = A3, EA4+CAMSIS = A4, Li D-comparable S1 = B1, EA4+Hill (Le) = B2, EA4+ISEI (Le) = B3); A1-EUR, A3-EUR, B1-pruned and MVP AFR were not computed", "outcome_precision_mvmr.tsv", rd("outcome_precision_mvmr.tsv")))
sheets.append(("joint_set_composition", "Composition of the 11 SES-versus-SES MVMR instrument sets (CAD Aragam EUR): SNPs from each exposure's own list and genome-wide significant per exposure", "mvmr_joint_set_composition.tsv", rd("mvmr_joint_set_composition.tsv")))
# Instrument lists: identifiers and alleles for every set. Effect sizes (beta, SE, as released) are included
# only for the SSGAC EA sets (EA4, EA3): SSGAC confirmed by email on 2026-10-11 that the effect allele, beta and
# SE of the instrument SNPs may be included. Other sets stay identifiers/alleles only until their terms are
# checked (Kweon is CC BY-NC-SA; others unconfirmed) - readers obtain those effect sizes from each source (S1).
BETA_OK = {"ea_okbay2022", "ea_lee2018"}
import glob
_inst = []
for f in sorted(glob.glob("data/instruments/*_instruments.tsv")):
    _ds = os.path.basename(f).replace("_instruments.tsv", "")
    d = pd.read_csv(f, sep="\t", usecols=["rsid", "chr", "pos", "effect_allele", "other_allele", "beta", "se"])
    if _ds not in BETA_OK:
        d["beta"] = pd.NA
        d["se"] = pd.NA
    d = d.rename(columns={"beta": "beta_as_released", "se": "se_as_released"})
    d.insert(0, "dataset", _ds)
    _inst.append(d)
_inst = pd.concat(_inst, ignore_index=True)
_inst.to_csv("manuscript/tables/instrument_lists.tsv", sep="\t", index=False)
sheets.append(("instrument_lists", "Instrument SNPs per exposure and mediator: rsID, GRCh37 position and alleles; beta and SE as released (per SD of educational attainment, effect allele as shown) for the SSGAC education sets (EA4, EA3) only, with SSGAC permission. For other sets effect sizes are available from each source GWAS (see S1)", "manuscript/tables/instrument_lists.tsv", _inst))
sheets.append(("gses_genetic_corr", "Cross-trait LDSC genetic correlation of the shared SES factor (gSES) with each SES measure (unconstrained estimates can exceed 1)", "gses_feasibility_rg.tsv", rd("gses_feasibility_rg.tsv")))
sheets = [(f"S{i}_{n}", t, s_, d) for i, (n, t, s_, d) in enumerate(sheets, 1)]
contents = pd.DataFrame([(n, t, "results/tables/" + s if not s.startswith(("data/", "manuscript/")) else s, len(d)) for n, t, s, d in sheets],
                        columns=["sheet", "title", "source file", "rows"])
with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    contents.to_excel(xw, sheet_name="Contents", index=False)
    for n, t, s, d in sheets:
        d.to_excel(xw, sheet_name=n[:31], index=False)
print(f"wrote {OUT}: {len(sheets)} tables")
