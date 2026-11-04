# van Zuydam 2021 PAD — parse-error diagnosis, row-count reconciliation, validation verdict

Performed 2026-09-18 on the files in `data/raw/`. Read-only: no file in `data/raw/` was
modified, moved, re-extracted or re-downloaded. All findings are reproducible from the
commands in the "How this was checked" section.

## 1. Diagnosis: the checks ran against the archive, not the summary statistics

The two `ERROR` strings in `results/logs/2026-06-30_01_integrity_summary.tsv` are a path
bug in `R/01_download_sumstats.R`, not file corruption.

The auto-download registry (line 68) names the target `pad_vanzuydam2021.tar.gz`. Phase B
builds the inventory with `path = target_path`, so Phase C verified the **tarball**. The run
log confirms it: `[verify] pad_vanzuydam2021  (pad_vanzuydam2021.tar.gz)`. The extraction
block (lines 243-253) unpacks the archive but never redirects the inventory at the inner file.

Each check then failed in its own characteristic way:

| Check | What happened |
| --- | --- |
| `count_rows()` | `path` matches `\\.gz$`, so it ran `gunzip -c pad_vanzuydam2021.tar.gz \| wc -l`, which decompresses the tar *stream* (tar headers plus seven still-compressed `.gz` members) and counts newline **bytes** in that binary blob: **7,349,497** lines, minus 1 for the assumed header = **7,349,496**. Reproduced exactly today. |
| `peek_columns()` | `data.table::fread()` on binary input: `ERROR: embedded nul in string: ...` |
| `pvalue_sanity()` | same `fread()` call on the same binary input: `ERROR: invalid multibyte string 1` |

No corruption anywhere: `gzip -t` passes on the archive **and** on all seven extracted members,
and every extracted member's byte size matches the `tar -tzvf` listing.

## 2. Row-count reconciliation

**10,250,115 is correct**; it is the data-row count of the extracted
`final_meta/pad_primary.txt.gz`. Verified two independent ways:

- `gunzip -c pad_primary.txt.gz | wc -l` = 10,250,116 lines, minus the header = 10,250,115
- `pandas.read_csv` parses the whole file to shape (10,250,115, 11)

**7,349,496 is not a measurement of any dataset** — it is the newline-byte count of a
compressed archive, as shown above.

The two numbers were never in conflict about the data. `data/raw/MANIFEST.tsv` already records
the correct values as committed in `75f6e59`: `n_snps = 10250115`, `sha256 = e4c9c17ecaac5808...`
(which I recomputed today on `pad_primary.txt.gz` and matched), the correct column order, and a
note preserving the tarball's own hash `f3b8c7b3...`. The project notes took their figure from the
manifest. The stale artifact is the integrity log alone.

## 3. Which file belongs in the primary analysis

**`data/raw/pad_vanzuydam2021/final_meta/pad_primary.txt.gz`** — the unstratified PAD
meta-analysis. The bundled data dictionary (`data-dictionary-GoLEAD-2020.txt`, which says
"6 files" and then lists 7) describes the others as diabetes and smoking strata plus two
interaction files carrying per-stratum columns and `p_int`. Those are subgroup material at most,
never the 4x4 grid cell.

| File | Data rows | Role |
| --- | --- | --- |
| `pad_primary.txt.gz` | 10,250,115 | **PRIMARY grid outcome** |
| `pad_diabetesvsnodiabetes.txt.gz` | 9,550,199 | interaction (`p_int`, per-stratum columns) |
| `pad_eversmokervsneversmoker.txt.gz` | 9,398,728 | interaction (`p_int`, per-stratum columns) |
| `pad_eversmoker.txt.gz` | 9,348,913 | stratum |
| `pad_diabetes.txt.gz` | 8,106,662 | stratum |
| `pad_nodiabetes.txt.gz` | 8,025,221 | stratum |
| `pad_neversmoker.txt.gz` | 7,426,645 | stratum |

Per-file sizes and sha256 in `pad_member_inventory.tsv`.

## 4. Validation verdict for `pad_primary.txt.gz`: passes

- **Structure:** 11 columns, matching the manifest note exactly; zero missing values in any column.
- **p-values:** all 10,250,115 parse and all fall in (0, 1]. Median 0.49, min 5e-17,
  q0.001 = 5.5e-4, 109 variants at p < 5e-8, 590 at p < 1e-5. A null-dominated distribution with a
  credible tail for a 12,086-case GWAS.
- **Positive control:** the eight strongest signals are all 9p21.3 (`rs4977575`, `rs1333046`,
  `rs10757278`, `rs1333047`, ...; OR 1.13, p 5e-17), the canonical atherosclerosis locus. The file
  is what it claims to be.
- **Sample size fields:** `ncases` max 12,086 (matches the published figure), `ncontrols` max
  449,553 against the manifest's 449,548. Per-variant N varies by design (5,010 distinct `ncases`
  values); the five-control difference is a per-variant accounting artefact, not a different study.
- **Effect/SE/EAF fields:** no non-positive `or` or `or_se`; EAF within [0, 1] with 107 variants at
  exactly 0 and 5 at exactly 1 (monomorphic in the meta — excludable at instrument selection).
- **Identifiers:** all 10,250,115 `snp` values unique; 92.5% rsIDs, the rest `chr:pos:I/D` indel
  identifiers. 700,853 duplicated `(chr, bp)` pairs, as expected for multi-allelic sites and indels.
- **Genome build:** positions are GRCh37/hg19 — 407 chr1 variants sit beyond the GRCh38 chr1
  length (max 249,229,799 vs GRCh38's 248,956,422), which rules GRCh38 out.

So acquisition **and** validation now both pass for this dataset. What had never been validated was
validated here; the earlier log entry measured the wrong file.

## 5. Three caveats that carry into Phase 2/3

**(a) The header order contradicts the bundled data dictionary.** The dictionary states column 5 =
`snp`, column 6 = `eaf`. `pad_primary.txt.gz` is actually
`chr bp or_allele nonor_allele eaf snp or or_se p ncases ncontrols` — `eaf` and `snp` swapped. All
six other members follow the dictionary. Any positional read built from the dictionary would
silently load rsIDs as allele frequencies. `MANIFEST.tsv` already has the true order; read by
column *name*, never position. (The two interaction files also misspell the column `nonor_alelle`.)

**(b) `or`, `or_se` and `eaf` are rounded to two decimal places in the deposit — `p` is not.**
The file carries 267 distinct `or` values, 87 distinct `or_se`, 101 distinct `eaf`; 38% of variants
have `or_se` exactly 0.02. Consequences for MR:

- ln(OR) takes only 267 distinct values, quantization step ~0.00995 near OR = 1.
- 13.2% of variants (1,348,741) have OR exactly 1.00, i.e. outcome beta exactly 0.
- Median |ln(OR)| across the file is 0.020 — about two quantization steps. Outcome-side betas for
  arbitrary exposure instruments are therefore quantized at roughly half their own magnitude.
- Relative quantization of `or_se` is ~17% at the median (0.005/se), ~25% at the 90th percentile.
- Reported `p` is inconsistent with ln(OR)/`or_se` because `p` comes from the unrounded
  meta-analysis: for the top hit, p = 5e-17 against 9.9e-10 implied; median discrepancy 1.6 log10
  units among p < 1e-5 variants. **`p` is the high-precision field; effect and SE are coarse.**

This does not invalidate the file, and it is a property of the deposit rather than of the download.
It does mean a decision is owed before PAD estimates are reported: use the rounded `or`/`or_se` as
deposited, or reconstruct the outcome beta and SE from `p` plus the sign of ln(OR). Either choice
belongs in `DEVIATIONS.md` and in the limitations, alongside the existing note that PAD is the
power-limiting cell.

**(c) Chromosome coding is mixed at the edges:** 1-22 plus `'23'` (13 rows), `'24'` (2 rows) and
`'X'` (3 rows) — 18 rows in total, negligible in volume but needing an explicit rule rather than a
silent coercion.

## 6. Recommended fixes (not applied)

1. `R/01_download_sumstats.R` — after the extraction block, override the inventory path for
   `pad_vanzuydam2021` to `pad_vanzuydam2021/final_meta/pad_primary.txt.gz` so the three checks
   verify the summary statistics, keeping the tarball sha256 as download provenance. **As written,
   re-running the script would write the archive-level numbers back into both the log and the
   manifest.**
2. `results/logs/2026-06-30_01_integrity_summary.tsv` — the row for `pad_vanzuydam2021` is the
   only stale record; supersede it rather than silently editing a dated run log.
3. Project notes — the 10,250,115 figure is correct; no change needed.
4. Log the rounding finding (b) and the chosen handling of it as a deviation.

Note for completeness: the same error scan found one further anomaly unrelated to PAD —
`tdi_panukb` returned `pval_check = no-pvalue-col`, because Pan-UKB ships `neglog10_pval_meta`
rather than a `p` column. That is a real gap in that file's validation, though a benign one.

## How this was checked

```
gzip -t pad_vanzuydam2021.tar.gz ; gzip -t final_meta/*.txt.gz      # integrity
tar -tzvf pad_vanzuydam2021.tar.gz                                  # member listing vs extracted sizes
gunzip -c final_meta/pad_primary.txt.gz | wc -l                     # 10,250,116
gunzip -c pad_vanzuydam2021.tar.gz     | wc -l                      # 7,349,497 (reproduces the log)
pandas.read_csv(pad_primary.txt.gz, sep='\t')                       # shape, NAs, ranges, quantiles
hashlib.sha256                                                      # per-member hashes
```
