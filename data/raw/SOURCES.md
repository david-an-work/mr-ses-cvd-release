# Sources

> One section per dataset. Updated as part of T1.1 sumstat reconnaissance (2026-05-07).
> When a file is downloaded for the first time, add: file path, sha256, row count, and
> any column-renaming notes. The MANIFEST.tsv row stays the source of truth for
> URL/access metadata; this file is the place for the *narrative* about each dataset.

This file preserves the May reconnaissance notes and later corrections. Some dataset roles, proposed analyses,
overlap statements and access assessments below were superseded as the study developed. For the current study,
consult `ANALYSIS_PLAN.md`, the manuscript and `data/raw/MANIFEST.tsv`; see `LICENSE-DATA.md` for the
release's redistribution and licensing scope.

---

## Exposures

### EA — Lee et al. 2018 (SSGAC EA3)
- **Citation:** Lee JJ, Wedow R, Okbay A, et al. *Nat Genet* 50:1112–1121 (2018).
- **What's released:** Excl-23andMe meta, N = 766,345 EUR. The full meta with 23andMe (N ≈ 1,131,881) is restricted by 23andMe to the lead-SNP table only.
- **URL:** https://www.thessgac.org/data — file `GWAS_EA_excl23andMe.txt` after click-through ToS.
- **Beta scale (corrected 2026-10-03):** SD units of educational attainment per allele, not years of schooling. Evidence: SE·√(2pqN) median 1.05 across instruments, and median instrument |β| 0.013 (≈2.7 weeks per allele at SD ≈ 4 years; ≈0.7 weeks if read as years). The published "weeks per allele" figures are SD-scale betas multiplied by the EA SD.
- **Sample-overlap caveat:** UKB ≈ 442k. Any UKB-derived CVD outcome inherits substantial overlap → MR-RAPS as sensitivity per `ANALYSIS_PLAN.md` §8.
- **Ancestry:** EUR only.
- **Within-family caveat:** sibship MR studies report direct EA effects ~40% smaller than population estimates, consistent with substantial assortative-mating + dynastic contributions. Discuss in Limitations.

### EA — Okbay et al. 2022 (SSGAC EA4)
- **Citation:** Okbay A, Wedow R, Kong A, et al. *Nat Genet* 54:437–449 (2022).
- **What's released:** Excl-23andMe meta, N = 765,283 EUR. The full N ≈ 3.04M meta (including 23andMe) is **restricted** — only ~10,000 lead SNPs are deposited publicly. Confirmed via SSGAC EA4 FAQ.
- **URL:** https://www.thessgac.org/data — file `EA4_additive_excl_23andMe.txt.gz`. Mirror: OSF project https://osf.io/uegqv/ .
- **Beta scale (corrected 2026-10-03):** SD units of educational attainment per allele (SE·√(2pqN) median 1.07), consistent with Lee 2018.
- **Q2 (open question) resolved:** access is **mixed**; the excl-23andMe file is the standard MR primary going forward. Project decision: use EA4 excl-23andMe as primary; Lee 2018 as a secondary/sensitivity comparator.
- **Direct vs population effects:** Okbay reports direct effects ≈ 56% of population effects.

### Income — Hill et al. 2019 (SENSITIVITY comparator from 2026-10-03; was primary)
- **Role change 2026-10-03:** demoted to sensitivity comparator when Kweon 2025 became primary (see below and `DEVIATIONS.md`).
- **Citation:** Hill WD et al. *Nat Commun* 10:5741 (2019).
- **URL:** https://datashare.ed.ac.uk/handle/10283/3441 (DOI 10.7488/ds/2643). Two files: `Household_Income_UKBiobank.txt.gz` (392 MB; UKB-only N=286,301) and `MTAG_household_Income.txt.gz` (151 MB; effective N≈397K, MTAG-augmented with EA + cognition).
- **Use:** UKB-only file as the clean instrument; the MTAG file mixes EA/cognition signal and is harder to interpret in MVMR.
- **Beta scale (corrected 2026-10-08):** NOT standardized log-income. Hill 2019 analysed residuals of the 1–5 household-income band scale (residualised on age, sex, 40 PCs, array, batch; GWAS in BGENIE) with no stated standardisation; the released SEs imply a scale of ≈0.61 (median SE·√(2pqN)), not 1 SD. Effects are reported per released unit.
- **Sign:** higher beta = higher income (no flip).
- **Ancestry:** EUR (UKB white British).

### Income — Kweon et al. 2025 (PRIMARY from 2026-10-03)
- **Citation:** Kweon H, Burik CAP, Ning Y, … Koellinger PD. *Nat Hum Behav* 9:794–805 (2025). DOI 10.1038/s41562-024-02080-7.
- **What it is:** "Income Factor": MTAG meta-analysis across four income measures (individual, occupational, household, parental). Effective N = 668,288, 1KG-EUR-like. 162 loci. rg = 0.92 (s.e. 0.006) with EA.
- **URL (corrected 2026-10-03):** the paper's data-availability DOI 10.62891/aac85602 resolves to DeSci Nodes dPID 149 (https://nodes.desci.com/dpid/149), which lists `DATA/INC_FACTOR_GWAS.txt` ("Income GWAS meta-analysis across measures") plus household / individual / occupational income files. Files are served from `ipfs.desci.com`. The 2026-05-07 recon looked only at SSGAC and OSF `rg8sh` and wrongly concluded the file was not openly downloadable.
- **Licence:** CC BY-NC-SA (data files). Anything derived from it and deposited (e.g. on Zenodo) must carry the same licence.
- **UKB overlap:** UK Biobank is one of the contributing cohorts, so overlap with all four primary outcome GWAS remains.
- **Interpretation:** includes parental income, so the factor carries more dynastic / indirect genetic signal than "own income"; the paper reports roughly one quarter of the PGI's population effect is direct.
- **File (downloaded and validated 2026-10-03):** `data/raw/income_kweon2025_incfactor.txt`, uncompressed TSV, 9,131,507 rows (matches the paper's Table 1), sha256 `f56992516f07b716dce6fb07bf1527a50aa7c0855ba76404e958673f18d3656b`. Columns `SNP CHR BP A1 A2 EAF BETA SE P N`; autosomes only; **GRCh37**.
- **Effect allele:** `A1` (the deposit's FAQ does not define columns, so this was confirmed empirically: 100% sign agreement with Hill 2019 across 3,176 SNPs genome-wide significant in both; `EAF` correlates 0.998 with the 1000G EUR frequency of `A1`).
- **Beta scale:** SD units of the Income Factor (SE·√(2pqN) has median 1.00). Per-SNP N 485,056–692,936 (median 655,729).
- **Status:** primary income exposure (user decision D1, 2026-10-03).
- **Contact:** correspondence to Abdellaoui, Hill or Koellinger (per paper).

### TDI — Pan-UKB phenotype 189-irnt
- **Citation:** Pan-UKBB consortium / Karjalainen, Martin, Finucane et al. (Broad/Pan-UKBB 2022).
- **URL:** https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-189-both_sexes-irnt.tsv.bgz (1.81 GB; tabix-indexed). Manifest: `phenotype_manifest.tsv.bgz` in same bucket.
- **Ancestry:** N = 439,198 total — EUR 420,035; CSA 8,864; AFR 6,624; EAS 2,698; AMR 977; MID failed QC.
- **Beta scale:** inverse-rank-normal-transformed Townsend score (SD units).
- **CRITICAL sign convention:** higher TDI = **more** deprivation. To match the project's "higher = higher SES" convention, multiply beta and signed Z by −1 in `02_select_instruments.R`. Document in `EXPOSURE_DIRECTIONS.tsv`.

### Occupational status — Akimova et al. 2024 (CANDIDATE — pending decision)
- **Citation:** Akimova ES, Breen R, Brazel DM, Mills MC. *Nat Hum Behav* (2024).
- **What it is:** GWAS of three occupational status scales (CAMSIS, SIOPS, ISEI) in UKB. N = 273,157 EUR.
- **URL / accessions:** GWAS Catalog GCST90446160 (CAMSIS), GCST90446162 (SIOPS), GCST90446163 (ISEI). Open.
- **If adopted:** adds an occupational SES dimension distinct from EA/income and TDI, expanding the exposure dimension to 5. Adds rows to the primary results grid (5 × 4 = 20 cells) and would require Bonferroni at α = 0.05/20.
- **Pending user decision** before adding.

---

## Outcomes

### CAD — Aragam et al. 2022
- **Citation:** Aragam KG, Jiang T, Goel A, et al. *Nat Genet* 54:1803–1815 (2022). PMID 36474045.
- **URLs:**
  - Multi-ancestry: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132315/GCST90132315_buildGRCh37.tsv (GCST90132315). N cases 210,842 / N controls 1,167,328.
  - EUR-only stratum: same FTP path with GCST90132314. N cases 181,522 / N controls 984,168.
- **Status:** still the most current openly downloadable multi-ancestry CAD GWAS as of 2026-05-07. Tcheandjieu MVP 2022 (PMID 35915156) is larger but only top-hits sumstats are deposited.
- **Ancestry:** ~85% EUR + ~14% EAS (BBJ) + small NR.
- **UKB inclusion:** yes (CARDIoGRAMplusC4D 1.0M meta includes UKB).
- **File size:** ~3.0 GB per file uncompressed; uploaded 2025-02.

### PAD — van Zuydam 2021 (PRIMARY) + FinnGen R12 (sensitivity)
- **Original plan** assumed Klarin et al. 2019 (Nat Med, MVP) at figshare 7811639. **Both elements were wrong:** figshare 7811639 hosts a van Zuydam CAD/T2D file, **not** Klarin PAD; Klarin 2019 PAD sumstats live on dbGaP `phs001672` (controlled-access, requires application) — incompatible with the project's "no biobank application" constraint.
- **Decision applied 2026-05-07: van Zuydam et al. 2021** *Circ Genom Precis Med* as primary PAD. Open figshare 12723404. EUR-only, 12,086 cases / 449,548 controls. Cleanest open alternative; smaller N than Klarin but citable without caveats.
- **FinnGen R12 PAD comparators** (Phase 5 sample-overlap sensitivity, public Google Cloud Storage URLs):
  - `I9_ATHSCLE` "Atherosclerosis, excluding cerebral and coronary sclerosis" — 19,807 cases / 463,106 controls. **Largest FinnGen PAD-adjacent endpoint; preferred FinnGen primary comparator.** URL: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_ATHSCLE.gz`
  - `I9_OTHPER` "Other peripheral vascular diseases" — 3,153 cases. Secondary FinnGen comparator. URL: `…/finngen_R12_I9_OTHPER.gz`
  - `I9_ARTEMBTHRLOW` "Arterial embolism / thrombosis lower extremity artery" — 1,297 cases. Tertiary; smallest but most specific. URL: `…/finngen_R12_I9_ARTEMBTHRLOW.gz`
  - All three: Finnish EUR; non-overlapping with UKB (clean from sample-overlap perspective with EA exposure). Manifest: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv`.
- **MVP-PheWAS Verma 2024 — REINSTATED 2026-10-03 as PAD replication (EUR + AFR; AMR descriptive).** Dropped on 2026-05-07 on the belief that dbGaP `phs002453` needed an eRA Commons account; that was wrong. `phs002453` is MVP's open-access accession ("publicly available results", no application); only `phs001672` (sensitive phenotypes / population subsets, incl. Klarin 2019 PAD) is controlled. Files: https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs002453/analyses/ (GIA and HARE ancestry schemes). Chosen trait: `Phe_440_2` "Atherosclerosis of the extremities" (GIA cases: all-ancestry 20,481; AFR 4,282; AMR 1,163; EUR 15,036; EAS not run). Files `MVP_R4.1000G_AGR.Phe_440_2.{META,AFR,AMR,EUR}.GIA.dbGaP.txt.gz` sit inside `phs002453.MVP_R4.1000G_AGR.GIA.PheCodes_CirculatorySystem_batch2.analysis-PI.MULTI.tar` (139 GB, not expanded on the server; HTTP range requests supported). EHR phecode definition, not a curated PAD phenotype. MVP does not include UKB. `Phe_443` "Peripheral vascular disease" (56,156 cases) is broader and mostly "unspecified"; not used. **Downloaded and validated 2026-10-03** as `data/raw/pad_mvp_phe440_2_{EUR,AFR,AMR,META}.txt.gz` (byte ranges of the tar; see `R/01_download_sumstats.R`). Columns `SNP_ID chrom pos ref alt ea af num_samples case_af num_cases control_af num_controls or ci pval r2 q_pval i2 direction`; `ea` = `alt`; **GRCh38** (rs10757278 at chr9:22,124,478), unlike the GRCh37 exposure files, so harmonize on rsID or lift over. **No SE column:** SE = (ln upper − ln lower)/3.92 from `ci`; reported p agrees with OR/CI. EUR: 19,709,222 rows, 1,569 at p<5e-8, top signals 9p21 (OR 0.84 for rs10757278-A, same direction as van Zuydam). AFR: 36,970,753 rows, none at p<5e-8 (top: chr15, and ABO at chr9:133 Mb). AMR: 27,867,029 rows, none at p<5e-8. META: 27,515,940 rows, 826 at p<5e-8. Checksums in `MANIFEST.tsv`. *Superseded text (2026-05-07):* Original plan would have used MVP-PheWAS (dbGaP `phs002453`, Verma et al. *Science* 2024) for ancestry-stratified PAD replication. Dropped because dbGaP "summary view" access requires an NIH eRA Commons account, which can only be created by an institutional Signing Official at Harvard's Office of Sponsored Programs (3-7 day bureaucratic loop). Project decision: keep the "no application" constraint and document the gap as an honest manuscript limitation: *"ancestry-stratified PAD replication was not feasible because the only large openly-described multi-ancestry PAD GWAS (MVP) is gated by NIH eRA Commons."*
- **No multi-ancestry post-Klarin PAD GWAS with openly downloadable full sumstats exists** as of 2026-05-07 (Q1 resolved negative).

### Heart failure — HERMES 2.0 (Henry/Levin et al. 2025) — PRIMARY
- **Citation:** Henry A, Levin MG, et al. (HERMES) *Nat Genet* (2025). PMID 40038546.
- **URL:** https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90726001-GCST90727000/GCST90726617/GCST90726617.tsv.gz (495 MB). Subtypes: GCST90726618 (ni-HF), 90726619 (ni-HFpEF), 90726620 (ni-HFrEF). Mirror set 90728695–698.
- **Sample size:** all-cause HF 153,174 cases / 1,793,175 controls. Subtypes: ni-HF 44,012; ni-HFrEF 5,406; ni-HFpEF 3,841 cases.
- **Ancestry:** EUR 91% / EAS 6.2% / AFR 2.2% / SAS 0.5% / AMR 0.1%.
- **Q3 (open question) resolved:** HERMES 2.0 supersedes Shah 2020 (47K cases) and Joseph 2022 (115K cases). Use HERMES 2.0 all-cause HF as primary; subtypes as supplementary.
- **Uploaded:** 2026-02. UKB included.

### Stroke — GIGASTROKE / Mishra et al. 2022
- **Citation:** Mishra A, Malik R, Hachiya T, et al. *Nature* 611:115–123 (2022). PMID 36180795.
- **URLs:**
  - Any-stroke multi: GCST90104534 — N cases 110,182 / N controls 1,503,898
  - Ischemic-stroke multi: GCST90104535
  - Any-stroke EUR: GCST90104539 — N cases 73,652 / N controls 1,234,808
  - Subtype + ancestry strata: GCST90104537 (SVS multi), 90104538 (LAS multi), 90104541–543 (CES/LAS/SVS EUR), 90104546–548 (EAS), 90104552–553 (AFR), 90104556–558 (HIS), 90104561–563 (SAS).
- **Ancestry:** EUR ~67% / EAS ~22% / AFR/HIS/SAS combined ~11%.
- **Avoid:** MTAG-augmented EUR strata (90162543–546). Inappropriate for MR — induced LD distortions and partial p-value sets.
- **Status:** no 2023–2026 stroke GWAS supersedes GIGASTROKE.

---

## Mediators

### BMI — Yengo et al. 2018 (GIANT + UKB)
- **URL:** https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files (also Zenodo record 1251813).
- **N:** ~700K EUR. Beta on IRNT-BMI (SD units).
- **Status:** still the practical default. The 2024 GIANT meta with 536 loci exists, but per-variant sumstats are not yet broadly released.

### Smoking — Wootton et al. 2020 (lifetime smoking index)
- **URL:** https://doi.org/10.5523/bris.10i96zb8gm0j81yz0q6ztei23d (Bristol data.bris).
- **N:** 462,690 EUR (UKB).
- **Beta scale (corrected 2026-10-08):** NOT SD units. The data.bris README states the betas are unstandardised (divide beta and SE by SD 0.6940093 for SD units); the released SEs imply ≈0.67. Smoking failed the mediation strength gate, so no smoking estimate is interpreted; indirect effects (a×b) are invariant to mediator scale.
- **Alternative:** Saunders et al. 2022 GSCAN multi-ancestry (UMN: https://conservancy.umn.edu/handle/11299/241912; ~3.4M total). Useful if ancestry-stratified mediation through smoking is desired; otherwise Wootton stays primary.

### Lipids — Graham et al. 2021 / GLGC
- **URL:** https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/trans_ancestry/
- Three files: LDL-C, HDL-C, triglycerides.
- **N:** 1.65M trans-ancestry (EUR 1.32M / EAS 147K / AFR 99K / HIS 48K / SAS 41K).
- **Units:** SD units — GLGC 2021 analysed inverse-normal-transformed lipid residuals, so betas are per SD, not mg/dL (corrected 2026-10-03; the earlier "mg/dL; TG log(mg/dL)" was wrong).

### SBP — Keaton et al. 2024 — PRIMARY (replaces Evangelou 2018)
- **Citation:** Keaton JM et al. *Nat Genet* (2024). PMID 38689001.
- **URL:** https://www.ebi.ac.uk/gwas/publications/38689001 — main accession GCST90310294 plus related.
- **N:** ~1,030,000 EUR (UKB + ICBP + MVP + BioVU).
- **Decision:** swap from Evangelou 2018 (757K) to Keaton 2024 (1.03M) — more recent, larger, same trait. Logged in DEVIATIONS.md as "newer release of same exposure" — not a methodological change.

### T2D — DIAMANTE / Mahajan et al. 2022
- **Citation:** Mahajan A et al. *Nat Genet* (2022).
- **URL:** https://diagram-consortium.org/downloads.html (free request form). Per-ancestry files for EUR/EAS/SAS/AFR/HIS available.
- **N:** 180,834 cases / 1,159,055 controls. 5-group trans-ancestry.
- **Status:** no 2024–2026 DIAMANTE update with public sumstats has superseded Mahajan 2022. T2DGGI 2024 sumstats are released through the same DIAGRAM portal but largely a partial re-meta with overlapping cohorts.

---

## Ancestry-stratified replication — what's actually possible

Realistically (as of 2026-05-07; PAD rows updated 2026-10-03), the (exposure × outcome × ancestry) cube populated with **publicly downloadable** sumstats:

- **EUR**: every exposure × every outcome — full coverage.
- **EAS**: BBJ covers CAD, HF, stroke, T2D, BMI, lipids; Pan-UKB EAS covers SES exposures but with N ≈ 2.7K — very low power for instrument selection.
- **AFR**: Pan-UKB AFR (~6.6K) plus GBMI AFR cover CAD/HF/stroke and BMI/lipids. SES instruments in AFR rest entirely on Pan-UKB AFR — very low power. **PAD AFR:** MVP `phs002453` `Phe_440_2` AFR (4,282 cases), open access; adopted 2026-10-03 (the earlier "no openly available AFR PAD" statement was wrong).
- **HIS / AMR / SAS**: GLGC and Pan-UKB give mediator coverage; outcomes only via Pan-UKB phecodes; SES exposures only via Pan-UKB. **PAD AMR:** MVP `Phe_440_2` AMR (1,163 cases), reported descriptively only. No open SAS PAD sumstats identified.

**Q4 (open question, All of Us) resolved negative** *(2026-09-06 update: the user holds All of Us Controlled Tier access; All of Us is deferred to a follow-up paper — see `DEVIATIONS.md`)*: as of 2026-05-07 All of Us has no openly downloadable sumstats — the All-by-All tables sit inside the Researcher Workbench Controlled Tier and require a workbench application + institutional DUC. AoU is therefore unusable under the project's "no application" constraint.

The realistic ancestry plan: **EUR** (full primary), **EAS** (CVD-strong, SES-weak), **AFR** (CVD-moderate, SES-very-weak), **HIS / SAS** (mostly mediator-only). The "ancestry equity gap" framing (`ANALYSIS_PLAN.md` §10) is well-justified by what's actually estimable.

---

## Sources to cite for the methods reconnaissance itself
- SSGAC EA4 FAQ (verified 2026-05-07): describes the with/without 23andMe split.
- Pan-UKB phenotype manifest (verified 2026-05-07).
- GWAS Catalog FTP browse for HERMES 2.0 (GCST907266** range, uploaded 2026-02).
- GWAS Catalog publication page for Keaton 2024 (PMID 38689001).
- All of Us Research Workbench documentation (verified 2026-05-07): All-by-All tables remain Controlled-Tier-only.

## Reference data (added 2026-10-03; stored in `data/reference/`, gitignored)

- **1000 Genomes LD reference panel (MRC-IEU, `1kg.v3.tgz`)** — PLINK bfiles for EUR (503 samples, 8,550,156 variants), AFR, AMR, EAS, SAS. Source: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz (downloaded 2026-10-03; 1,567,805,147 bytes; sha256 `06125160770fc65adf0296fe9a461011a910338bd106b4e319add13e81ab3670`). Extracted to `data/reference/1kg_v3/`. Used for local LD clumping (`ieugwasr::ld_clump(bfile=, plink_bin=)`), since the OpenGWAS API is not reachable from the analysis sandbox.
- **PLINK binary** — PLINK v1.90b6.10 (17 Jun 2019), bundled in the `plinkbinr` R package (GitHub `explodecomputer/plinkbinr`, recorded in `renv.lock`).
- **European LD scores (`eur_w_ld_chr`)** — 1000G EUR HapMap3 LD scores and weights from the LDSC developers, plus `w_hm3.snplist`. The original Broad URL (data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2) now redirects to a 404, so the copy at Zenodo record 8182036 (doi:10.5281/zenodo.8182036, "a gzipped copy of the European LD scores … provided by the Alkes Group") was used. Downloaded 2026-10-03; md5 `e2f16343c4cfaa76caa7d0c03d26b489` (matches the Zenodo record); sha256 `9537f00eb0d163a935aaa2cf04b358b7cf21852279b9c7925802526f6060b069`. Extracted to `data/reference/eur_w_ld_chr/`. Used for LD score regression genetic correlations (GenomicSEM `ldsc()`; T2.11).
