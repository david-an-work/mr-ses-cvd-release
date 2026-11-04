# Integrity-check addendum — 2026-10-03

Supersedes specific rows of `2026-06-30_01_integrity_summary.tsv`, which is left unedited as the
dated record of that run, and records the checks for the Kweon and MVP files downloaded on 2026-10-03. See `DEVIATIONS.md` (2026-10-03, T1.3).

| dataset | file checked | n_rows | sha256 | columns | p-value check | status |
| --- | --- | --- | --- | --- | --- | --- |
| pad_vanzuydam2021 | `pad_vanzuydam2021/final_meta/pad_primary.txt.gz` | 10,250,115 | e4c9c17ecaac5808da90ee91f726328359d85e234c50ef7b13ac027d9c00fd08 | chr, bp, or_allele, nonor_allele, eaf, snp, or, or_se, p, ncases, ncontrols | all 10,250,115 p in (0,1]; median 0.49; min 5e-17; 109 at p<5e-8; top signals 9p21.3 | PASS (validated 2026-09-18) |
| tdi_panukb | `tdi_panukb_189_irnt.tsv.bgz` (`neglog10_pval_meta`) | 28,987,534 | unchanged | as logged 2026-06-30 | 428,158 NA; p median 0.475; min 2.7e-13; 71 at p<5e-8 (10 chromosomes); 2,738 at p<1e-5; reported −log10 p identical to that implied by beta_meta/se_meta; λGC 1.12 | PASS (checked 2026-10-03) |
| income_kweon2025 | `income_kweon2025_incfactor.txt` | 9,131,507 | f56992516f07b716dce6fb07bf1527a50aa7c0855ba76404e958673f18d3656b | SNP, CHR, BP, A1, A2, EAF, BETA, SE, P, N | all P in (0,1]; median 0.416; min 3.2e-37; 13,283 at p<5e-8; P = f(BETA/SE) | PASS (downloaded + checked 2026-10-03) |
| pad_mvp_phe440_2_eur | `pad_mvp_phe440_2_EUR.txt.gz` | 19,709,222 | b0f251e89d3ddd98b9dae27d10065a3dcc0cdece0c2e43463f33cb7ab5cc6782 | SNP_ID, chrom, pos, ref, alt, ea, af, …, or, ci, pval, … | no NA; p median 0.468; min 5e-54; 1569 at p<5e-8; p consistent with OR/CI | PASS (downloaded + checked 2026-10-03) |
| pad_mvp_phe440_2_afr | `pad_mvp_phe440_2_AFR.txt.gz` | 36,970,753 | 7ec8c0fc813f20c16acb7f32b87f4e688c50db14bcd1186c1d32d3ee8746f591 | SNP_ID, chrom, pos, ref, alt, ea, af, …, or, ci, pval, … | no NA; p median 0.493; min 5.4e-08; 0 at p<5e-8; p consistent with OR/CI | PASS (downloaded + checked 2026-10-03) |
| pad_mvp_phe440_2_amr | `pad_mvp_phe440_2_AMR.txt.gz` | 27,867,029 | cb2128850220edaeb74d68df0f4f28131aad5dac637094637877d2c770ee6c9c | SNP_ID, chrom, pos, ref, alt, ea, af, …, or, ci, pval, … | no NA; p median 0.486; min 1.4e-07; 0 at p<5e-8; p consistent with OR/CI | PASS (downloaded + checked 2026-10-03) |
| pad_mvp_phe440_2_meta | `pad_mvp_phe440_2_META.txt.gz` | 27,515,940 | 4060e3e6ee4546133d8ef6cbc16864a9f63ae3c11a8d008492a6d9719418bff5 | SNP_ID, chrom, pos, ref, alt, ea, af, …, or, ci, pval, … | no NA; p median 0.496; min 5.7e-44; 826 at p<5e-8; p consistent with OR/CI | PASS (downloaded + checked 2026-10-03) |

Notes
- Kweon effect allele confirmed as `A1` and beta scale as SD units (see `SOURCES.md`). MVP files are GRCh38
  and give OR + 95% CI with no SE column; exposure files are GRCh37.
- The first download attempt was interrupted by dropped connections; files were completed by
  appending byte ranges and checked against their exact expected sizes (and gzip CRC for MVP).
- 2026-06-30 values for `pad_vanzuydam2021` (`n_rows = 7349496`; column/p-value `ERROR`) were
  measured on the `.tar.gz` archive, not the summary statistics. No corruption: `gzip -t` passes on
  the archive and all seven members. Archive sha256 (download provenance):
  f3b8c7b3373e01dee4c06e286e99c950f6251d153b38418d3c22060abac3be1c.
- The TDI count of genome-wide significant variants is pre-clumping; independent-instrument counts
  come from Phase 2 (T2.6, T2.10).
- Full PAD diagnosis: `results/logs/2026-09-18_PAD_VALIDATION.md`.
