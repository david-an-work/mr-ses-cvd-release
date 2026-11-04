"""03e_verify_proxies.py

Independent check of the proxy sensitivity inputs written by R/03d (separate language; no PLINK; shares no code with R).
  1. Reads the 1000 Genomes EUR (and AFR) PLINK .bed genotypes directly and, for every target/proxy pair in
     results/tables/proxy_log.tsv, computes the genotype (allele-dosage) correlation. Checks: same chromosome,
     distance <= 500 kb, proxy not palindromic, genotype r2 >= 0.8 and within 1e-4 of the logged PLINK r2, and (MVP AFR)
     genotype r2 >= 0.8 in AFR with the same phase sign. The phase sign decides which target allele the proxy's outcome effect allele tags.
  2. Re-reads each proxy SNP's outcome row from the raw file with the R/03c parsers (rule-1 key, rule 3).
  3. Re-derives the alignment of every proxy-augmented cell (univariable_proxy and mvmr_proxy) from the raw rows
     (original rows for non-proxied SNPs; independently phase-mapped proxy rows otherwise) and checks mr_keep and
     beta_y_aligned exactly, plus beta/SE/EAF carried from the raw proxy row.
Exits with status 1 on any mismatch. No MR estimate is computed.
Usage: python3 R/03e_verify_proxies.py   (from the repository root)
"""
import glob, importlib.util, os, sys
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("v03c", "R/03c_verify_source_rows.py")
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
H = v.H; REF = "data/reference/1kg_v3"
plog = pd.read_csv("results/tables/proxy_log.tsv", sep="\t")


class Bed:
    def __init__(self, prefix, need):
        bim = pd.read_csv(prefix + ".bim", sep="\t", header=None, usecols=[0, 1, 3, 4, 5],
                          names=["chr", "rsid", "pos", "a1", "a2"], dtype={"chr": str, "a1": str, "a2": str})
        bim["idx"] = np.arange(len(bim))
        self.bim = bim[bim.rsid.isin(need)].drop_duplicates("rsid").set_index("rsid")
        n = sum(1 for _ in open(prefix + ".fam")); self.n = n; self.nb = (n + 3) // 4
        self.f = open(prefix + ".bed", "rb"); assert self.f.read(3) == b"\x6c\x1b\x01"

    def dosage(self, rsid):  # dosage of bim a1
        self.f.seek(3 + int(self.bim.at[rsid, "idx"]) * self.nb)
        b = np.frombuffer(self.f.read(self.nb), dtype=np.uint8)
        g = np.stack([(b >> s) & 3 for s in (0, 2, 4, 6)], axis=1).ravel()[: self.n]
        return np.select([g == 0, g == 2, g == 3], [2.0, 1.0, 0.0], np.nan)

    def corr(self, t, p, t_allele, p_allele):
        """Correlation between dosage of t_allele at t and dosage of p_allele at p."""
        x, y = self.dosage(t), self.dosage(p)
        if self.bim.at[t, "a1"] != t_allele: x = 2 - x
        if self.bim.at[p, "a1"] != p_allele: y = 2 - y
        ok = ~np.isnan(x) & ~np.isnan(y)
        return float(np.corrcoef(x[ok], y[ok])[0, 1])


need = set(plog.rsid) | set(plog.proxy)
eur = Bed(os.path.join(REF, "EUR"), need); afr = Bed(os.path.join(REF, "AFR"), need)
fail = []

# raw outcome rows for the proxy SNPs, keyed exactly as in R/03c (GRCh37 position from the EUR .bim)
pref = eur.bim.loc[sorted(set(plog.proxy))].reset_index()[["rsid", "chr", "pos", "a1", "a2"]]
pref["key_pair"] = [frozenset((a, b)) for a, b in zip(pref.a1, pref.a2)]
v.inst = pref.set_index("rsid")
mapped = {}  # outcome -> DataFrame of proxy rows mapped to the target SNP
checks = []
for o, g in plog.groupby("outcome"):
    raw_p, _ = v.raw_rows(o)
    rows = []
    for _, r in g.iterrows():
        t, p = r.rsid, r.proxy
        if p not in raw_p.index: fail.append((o, t, p, "proxy row not unique/ok in raw file")); continue
        y = raw_p.loc[p]
        tb, pb = eur.bim.loc[t], eur.bim.loc[p]
        c = eur.corr(t, p, tb.a1, y.ea)  # sign: does target a1 travel with the proxy's outcome effect allele?
        ea_t, oa_t = (tb.a1, tb.a2) if c > 0 else (tb.a2, tb.a1)
        rec = dict(outcome=o, rsid=t, proxy=p, dist=abs(int(tb.pos) - int(pb.pos)), same_chr=tb.chr == pb.chr,
                   proxy_pal=(pb.a1, pb.a2) in v.PAL, r2_geno_eur=c * c, r2_plink=r.proxy_r2, ea_t=ea_t)
        if o == "pad_mvp_afr":
            ca = afr.corr(t, p, tb.a1, y.ea) if (t in afr.bim.index and p in afr.bim.index) else np.nan
            rec.update(r2_geno_afr=ca * ca, afr_same_phase=bool(np.sign(ca) == np.sign(c)))
            if not rec["afr_same_phase"] or not ca * ca >= 0.8: fail.append((o, t, p, "AFR phase/r2"))
        if not (rec["same_chr"] and rec["dist"] <= 500_000 and not rec["proxy_pal"] and c * c >= 0.8 and abs(c * c - r.proxy_r2) < 1e-4):
            fail.append((o, t, p, "pair check"))
        checks.append(rec)
        rows.append(dict(rsid=t, ea=ea_t, oa=oa_t, eaf=y.eaf, beta=y.beta, se=y.se))
    mapped[o] = pd.DataFrame(rows).set_index("rsid")

# re-derive every proxy-augmented cell
spec2 = importlib.util.spec_from_file_location("v03c_orig", "R/03c_verify_source_rows.py")
vo = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(vo)
n_cells = 0; n_proxy_rows = 0
for o in sorted(set(os.path.basename(f).split("__")[1][:-7] for f in glob.glob(os.path.join(H, "univariable_proxy", "*.tsv.gz")))):
    raw, _ = vo.raw_rows(o)
    raw_all = pd.concat([raw[["ea", "oa", "eaf", "beta", "se"]], mapped.get(o, pd.DataFrame(columns=raw.columns))])
    assert raw_all.index.is_unique
    tol = 0.085 if o == "pad_vanzuydam" else 0.08
    for f in sorted(glob.glob(os.path.join(H, "univariable_proxy", f"*__{o}.tsv.gz")) + glob.glob(os.path.join(H, "mvmr_proxy", f"*__{o}.tsv.gz"))):
        c = pd.read_csv(f, sep="\t", dtype={"ea_x": str, "oa_x": str, "ea_y": str, "oa_y": str, "proxy": str})
        k_ok, b_ok, _ = v.align(c, raw_all, tol, o == "pad_mvp_afr")
        pr = c[c.proxy.notna()].set_index("rsid")
        if len(pr):
            m = mapped[o].loc[pr.index]
            vals = (pr.ea_y.str.upper() == m.ea).all() and v.close(pr.beta_y, m.beta) and v.close(pr.se_y, m.se) and v.close(pr.eaf_y, m.eaf)
            n_proxy_rows += len(pr)
        else:
            vals = True
        n_cells += 1
        if not (k_ok and b_ok and vals): fail.append((o, os.path.basename(f), k_ok, b_ok, vals))

ck = pd.DataFrame(checks)
ck.to_csv("results/tables/proxy_verify.tsv", sep="\t", index=False)
print("pairs checked:", len(ck), "| max |r2_geno - r2_plink|:", round((ck.r2_geno_eur - ck.r2_plink).abs().max(), 4),
      "| min genotype r2:", round(ck.r2_geno_eur.min(), 3), "| max distance:", ck.dist.max())
print("cells checked:", n_cells, "| proxied rows checked:", n_proxy_rows)
if fail:
    print("FAILURES:", fail[:20], len(fail)); sys.exit(1)
print("All proxy pairs, phase mappings and proxy-augmented cells re-derive independently.")
