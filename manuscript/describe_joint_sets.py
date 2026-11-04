"""describe_joint_sets.py

Descriptive composition of the 11 SES-versus-SES MVMR instrument sets (CAD Aragam EUR harmonized cells; no estimation).
For each model and exposure: the number of joint-set SNPs that are in that exposure's own univariable instrument list
(the joint set is the union of the component lists re-clumped at r2 < 0.001 on the minimum p across exposures, plan §6),
and the number genome-wide significant (p < 5e-8, from the harmonized betas/SEs) for that exposure; plus SNPs in two or
more own lists and SNPs genome-wide significant for two or more exposures.
Usage (repository root): python3 manuscript/describe_joint_sets.py  ->  results/tables/mvmr_joint_set_composition.tsv
"""
import pandas as pd
from scipy.stats import norm

H = "data/harmonized/large/mvmr"
INST = {"EA": "ea_okbay2022", "Income": "income_kweon2025", "TDI": "tdi_panukb", "Occupation": "occupation_akimova2024",
        "Income_Hill": "income_hill2019", "TDI_EUR": "tdi_panukb_eur", "Occupation_ISEI": "occupation_akimova2024_isei"}
own = {k: set(pd.read_csv(f"data/instruments/{v}_instruments.tsv", sep="\t").rsid) for k, v in INST.items()}
app = pd.read_csv("results/tables/mvmr_appraisal.tsv", sep="\t")
models = app.loc[app.outcome == "cad_aragam_eur", ["model", "model_label"]].drop_duplicates()
rows = []
for m, lab in models.itertuples(index=False):
    d = pd.read_csv(f"{H}/{m}__cad_aragam_eur.tsv.gz", sep="\t")
    d = d[d.status == "ok"] if "status" in d else d
    exps = app.loc[(app.model == m) & (app.outcome == "cad_aragam_eur"), "exposure"].tolist()
    def col(e, pref):
        for c in (f"{pref}_{INST[e]}", f"{pref}_{e}", f"{pref}_Income" if e == "Income_Hill" else None):
            if c and c in d.columns:
                return c
        raise KeyError((m, e))
    gw = pd.DataFrame({e: 2 * norm.sf(abs(d[col(e, "b")] / d[col(e, "se")])) < 5e-8 for e in exps})
    mem = pd.DataFrame({e: d.rsid.isin(own[e]).values for e in exps})
    for e in exps:
        rows.append({"model": m, "model_label": lab, "n_joint_snps": len(d), "exposure": e,
                     "own_list_size": len(own[e]), "in_own_list": int(mem[e].sum()), "gw_significant": int(gw[e].sum()),
                     "model_snps_in_2plus_own_lists": int((mem.sum(1) >= 2).sum()),
                     "model_snps_gw_for_2plus_exposures": int((gw.sum(1) >= 2).sum())})
R = pd.DataFrame(rows)
R.to_csv("results/tables/mvmr_joint_set_composition.tsv", sep="\t", index=False)
print(len(R))
