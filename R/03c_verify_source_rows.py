"""03c_verify_source_rows.py

Independent check of R/03 (separate language and parser; shares no code with R/03).
For each outcome file:
  1. Re-reads the raw file with its own column rules and finds, for every union instrument SNP, the rows that match
     the frozen rule-1 key (GRCh37 files: chr + GRCh37 position + unordered allele pair; MVP and FinnGen: rsID +
     unordered allele pair). Exactly one matching row = ok (rule 3).
  2. Compares the R/03 outcome extract with the raw rows: status, alleles, beta (ln OR where applicable), SE (MVP:
     from the 95% CI; van Zuydam: or_se / OR) and EAF.
  3. Re-derives the alignment of every harmonized cell (univariable and MVMR; mediation cells use the R/03 mediator
     extracts) from the raw values and the exposure alleles: sign from the allele mapping, palindrome flip and
     ambiguity from the frequency rule (tolerance 0.08; van Zuydam 0.085; MVP AFR drops all palindromes), and checks
     mr_keep and beta_y_aligned exactly.
Exits with status 1 on any mismatch. No MR estimate is computed.
Usage: python R/03c_verify_source_rows.py   (from the repository root)
"""
import glob, gzip, math, os, sys
import numpy as np
import pandas as pd

RAW = "data/raw"; H = "data/harmonized/large"
PAL = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

# union of instrument SNPs with GRCh37 positions (from the instrument files, as selected by R/02)
inst = pd.concat([pd.read_csv(f, sep="\t", usecols=["rsid", "chr", "pos", "effect_allele", "other_allele"], dtype={"chr": str})
                  for f in sorted(glob.glob("data/instruments/*_instruments.tsv"))]).drop_duplicates("rsid")
inst["key_pair"] = [frozenset((a.upper(), b.upper())) for a, b in zip(inst.effect_allele, inst.other_allele)]
inst = inst.set_index("rsid")

def ci_se(lo, hi):
    return (np.log(hi.astype(float)) - np.log(lo.astype(float))) / (2 * 1.959963984540054)

def raw_rows(o):
    """Return DataFrame: rsid (instrument), ea, oa, beta, se, eaf (one or more rows per rsid before rule 3)."""
    if o in ("cad_aragam_eur", "cad_aragam_multi", "stroke_gigastroke_multi", "stroke_gigastroke_eur", "hf_hermes2"):
        f = {"cad_aragam_eur": "cad_aragam2022_eur.tsv", "cad_aragam_multi": "cad_aragam2022_multi.tsv",
             "stroke_gigastroke_multi": "stroke_gigastroke_any_multi.tsv.gz", "stroke_gigastroke_eur": "stroke_gigastroke_any_eur.tsv.gz",
             "hf_hermes2": "hf_hermes2_2025_hfall.tsv.gz"}[o]
        d = pd.read_csv(os.path.join(RAW, f), sep="\t", usecols=["chromosome", "base_pair_location", "effect_allele", "other_allele",
                        "beta", "standard_error", "effect_allele_frequency"], dtype={"chromosome": str})
        d = d.rename(columns={"chromosome": "chr", "base_pair_location": "pos", "effect_allele": "ea", "other_allele": "oa",
                              "standard_error": "se", "effect_allele_frequency": "eaf"})
        return by_position(d)
    if o == "pad_vanzuydam":
        d = pd.read_csv(os.path.join(RAW, "pad_vanzuydam2021", "final_meta", "pad_primary.txt.gz"), sep="\t", dtype={"chr": str})
        d = d.rename(columns={"bp": "pos", "or_allele": "ea", "nonor_allele": "oa"})
        d["beta"] = np.log(d["or"].astype(float)); d["se"] = d["or_se"].astype(float) / d["or"].astype(float)
        return by_position(d[["chr", "pos", "ea", "oa", "eaf", "beta", "se"]])
    if o in ("cad_mvp_phe411_eur", "pad_mvp_eur", "pad_mvp_afr"):
        f = {"cad_mvp_phe411_eur": "cad_mvp_phe411_EUR.txt.gz", "pad_mvp_eur": "pad_mvp_phe440_2_EUR.txt.gz",
             "pad_mvp_afr": "pad_mvp_phe440_2_AFR.txt.gz"}[o]
        d = pd.read_csv(os.path.join(RAW, f), sep="\t", usecols=["SNP_ID", "ref", "alt", "ea", "af", "or", "ci"])
        d = d[d.SNP_ID.isin(inst.index)].copy()
        assert (d.ea == d.alt).all()
        lo, hi = d.ci.str.split(",", expand=True)[0], d.ci.str.split(",", expand=True)[1]
        out = pd.DataFrame({"rsid": d.SNP_ID, "ea": d.alt, "oa": d.ref, "eaf": d.af.astype(float),
                            "beta": np.log(d["or"].astype(float)), "se": ci_se(lo, hi)})
        return by_rsid(out)
    if o in ("hf_mvp_phe428_eur", "stroke_mvp_phe433_21_eur"):
        f = {"hf_mvp_phe428_eur": "hf_mvp_phe428_EUR.gwascat.tsv.gz", "stroke_mvp_phe433_21_eur": "stroke_mvp_phe433_21_EUR.gwascat.tsv.gz"}[o]
        d = pd.read_csv(os.path.join(RAW, f), sep="\t", usecols=["rsid", "effect_allele", "other_allele", "effect_allele_frequency",
                        "odds_ratio", "ci_lower", "ci_upper"], na_values=["#NA", "NA"])
        d = d[d.rsid.isin(inst.index)].copy()
        out = pd.DataFrame({"rsid": d.rsid, "ea": d.effect_allele, "oa": d.other_allele, "eaf": d.effect_allele_frequency.astype(float),
                            "beta": np.log(d.odds_ratio.astype(float)), "se": ci_se(d.ci_lower, d.ci_upper)})
        return by_rsid(out)
    if o == "pad_finngen_athscle":
        d = pd.read_csv(os.path.join(RAW, "pad_finngen_r12_athscle.tsv.gz"), sep="\t", usecols=["ref", "alt", "rsids", "beta", "sebeta", "af_alt"])
        d = d[d.rsids.notna()].copy()
        d["rsid"] = d.rsids.str.split(",")
        d = d.explode("rsid")
        d = d[d.rsid.isin(inst.index)]
        return by_rsid(pd.DataFrame({"rsid": d.rsid, "ea": d.alt, "oa": d.ref, "eaf": d.af_alt, "beta": d.beta, "se": d.sebeta}))
    raise ValueError(o)

def by_position(d):
    d = d.copy(); d["chr"] = d["chr"].astype(str).str.replace("chr", "", regex=False); d["pos"] = d["pos"].astype(np.int64)
    ref = inst.reset_index()[["rsid", "chr", "pos", "key_pair"]].copy(); ref["pos"] = ref["pos"].astype(np.int64)
    m = d.merge(ref, on=["chr", "pos"])
    return finish(m)

def by_rsid(d):
    m = d.merge(inst.reset_index()[["rsid", "key_pair"]], on="rsid")
    return finish(m)

def finish(m):
    m = m.copy(); m["ea"] = m.ea.astype(str).str.upper(); m["oa"] = m.oa.astype(str).str.upper()
    m = m[[frozenset((a, b)) == k for a, b, k in zip(m.ea, m.oa, m.key_pair)]]
    n = m.groupby("rsid").size()
    return m[m.rsid.isin(n[n == 1].index)].set_index("rsid")[["ea", "oa", "eaf", "beta", "se"]], set(n[n > 1].index)

def close(a, b, rtol=1e-9):
    a = np.asarray(a, float); b = np.asarray(b, float)
    both_nan = np.isnan(a) & np.isnan(b)
    return np.all(both_nan | (np.abs(a - b) <= rtol * np.maximum(1, np.abs(b))))

def align(cell, raw, tol, drop_all_pal):
    """Independent re-derivation of mr_keep and the aligned beta for one harmonized cell."""
    c = cell.set_index("rsid")
    exp_keep = pd.Series(False, index=c.index); exp_beta = pd.Series(np.nan, index=c.index)
    for r in c.index:
        if r not in raw.index:
            continue
        ea, oa, fx = str(c.at[r, "ea_x"]).upper(), str(c.at[r, "oa_x"]).upper(), c.at[r, "eaf_x"]
        y = raw.loc[r]
        if (y.ea, y.oa) == (ea, oa): b, fy = y.beta, y.eaf
        elif (y.ea, y.oa) == (oa, ea): b, fy = -y.beta, (1 - y.eaf if pd.notna(y.eaf) else np.nan)
        else: raise AssertionError(f"allele pair mismatch {r}")
        pal = (ea, oa) in PAL
        fa = 0.5 if pd.isna(fx) else fx; fb = 0.5 if pd.isna(fy) else fy
        if pal and ((fa < 0.5 and fb > 0.5) or (fa > 0.5 and fb < 0.5)): b = -b
        amb = pal and ((0.5 - tol < fa < 0.5 + tol) or (0.5 - tol < fb < 0.5 + tol))
        keep = not amb and not (drop_all_pal and pal)
        exp_keep[r] = keep; exp_beta[r] = b
    ok_keep = (exp_keep == c.mr_keep.astype(bool)).all()
    k = exp_keep[exp_keep].index
    ok_beta = close(c.loc[k, "beta_y_aligned"], exp_beta[k], rtol=1e-9)
    return ok_keep, ok_beta, int(exp_keep.sum())

def main():
    fail = []; rows = []
    outcomes = sorted(os.path.basename(f)[:-7] for f in glob.glob(os.path.join(H, "outcome_extracts", "*.tsv.gz")))
    for o in outcomes:
        raw, dups = raw_rows(o)
        e = pd.read_csv(os.path.join(H, "outcome_extracts", o + ".tsv.gz"), sep="\t").set_index("rsid")
        ok_r = e.index[e.status == "ok"]
        st_ok = set(ok_r) == set(raw.index) and set(e.index[e.status == "duplicate"]) == dups
        common = raw.index.intersection(ok_r)
        vals_ok = (e.loc[common, "ea_y"].str.upper() == raw.loc[common, "ea"]).all() and (e.loc[common, "oa_y"].str.upper() == raw.loc[common, "oa"]).all() \
            and close(e.loc[common, "beta_y"], raw.loc[common, "beta"]) and close(e.loc[common, "se_y"], raw.loc[common, "se"]) \
            and close(e.loc[common, "eaf_y"], raw.loc[common, "eaf"])
        tol = 0.085 if o == "pad_vanzuydam" else 0.08
        n_cells = 0
        for f in sorted(glob.glob(os.path.join(H, "univariable", f"*__{o}.tsv.gz")) + glob.glob(os.path.join(H, "mvmr", f"*__{o}.tsv.gz"))):
            c = pd.read_csv(f, sep="\t", dtype={"ea_x": str, "oa_x": str})
            k_ok, b_ok, n_keep = align(c, raw, tol, o == "pad_mvp_afr")
            n_cells += 1
            if not (k_ok and b_ok): fail.append((o, os.path.basename(f), k_ok, b_ok))
        rows.append(dict(outcome=o, union_ok_raw=len(raw), union_ok_extract=len(ok_r), duplicates=len(dups),
                         status_match=st_ok, values_match=bool(vals_ok), cells_checked=n_cells))
        if not (st_ok and vals_ok): fail.append((o, "extract", st_ok, vals_ok))
        print(o, len(raw), "ok rows;", n_cells, "cells;", "status", st_ok, "values", bool(vals_ok), flush=True)

    # mediation cells: alignment logic re-derived from the R/03 mediator extracts (mediator GWAS read by the registry readers)
    for f in sorted(glob.glob(os.path.join(H, "mediation", "*.tsv.gz"))):
        c = pd.read_csv(f, sep="\t", dtype={"ea_x": str, "oa_x": str, "ea_y": str, "oa_y": str})
        okc = c[c.status == "ok"]
        raw = pd.DataFrame({"ea": okc.ea_y.str.upper().values, "oa": okc.oa_y.str.upper().values, "eaf": okc.eaf_y.values,
                            "beta": okc.beta_y.values, "se": okc.se_y.values}, index=okc.rsid.values)
        k_ok, b_ok, _ = align(c, raw, 0.08, False)
        if not (k_ok and b_ok): fail.append(("mediation", os.path.basename(f), k_ok, b_ok))
    pd.DataFrame(rows).to_csv("results/tables/harmonization_verify_source.tsv", sep="\t", index=False)
    if fail:
        print("FAILURES:", fail); sys.exit(1)
    print("All outcome extracts match the raw files, and every harmonized cell's alignment re-derives exactly.")


if __name__ == "__main__":
    main()
