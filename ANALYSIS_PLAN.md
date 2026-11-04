# Analysis Plan — version 1.0 (frozen 2026-10-07)

> **Status: frozen.** Approved by the study author on 2026-10-07. This version is draft v0.6 as approved, after the
> whole-plan consistency review (`docs/prereg/plan_v0.6_review.md` v2). It replaces v0.1 (archived at
> `docs/archive/ANALYSIS_PLAN_v0.1.md`); every change from v0.1 is logged in `results/logs/DEVIATIONS.md` (entries dated
> 2026-10-07, plus earlier dataset entries). It is to be registered on OSF **before any Phase 4 causal estimate** is
> computed; the OSF registration identifier will be recorded in `README.md` and `TASKS.md` (T9.8), not by editing this
> file. Any later change is a dated deviation, logged before the results it affects are examined.
>
> **Data already examined before the freeze (no causal estimate has been computed for any cell):** exposure-side
> instrument selection and strength diagnostics (`R/02`–`R/02e`), genetic correlations and LDSC intercepts among
> exposures and mediators, a simulation (`R/02f`), `qhet_mvmr` checks on simulated data (`R/02h`), outcome-blind
> precision from outcome-GWAS standard errors only (`R/02i`), the gSES feasibility check (`R/02j`), and outcome-file
> integrity checks: allele and imputation audits at the instrument SNPs, an SE-scale check of the rounded van Zuydam
> file (`R/02g`), and a value-by-value identity check of two copies of MVP `Phe_411`. Some of these checks read
> SNP-level outcome values (van Zuydam OR and p at selected variants; `Phe_411` OR/CI at instrument SNPs, compared only
> for identity between copies). No SNP-level outcome associations were combined into any MR estimate.
>
> **Supporting documents** (archived unchanged in `docs/prereg/`, checksums in `osf_registration/prereg_manifest.tsv`):
> `HARMONIZATION_RULES_draft.md` draft v3; `D7_overlap_comparison.md` draft v4 with `d7_evidence.tsv`;
> `reconstruction_feasibility.md` v2; `pilot_lid_cad_spec.md` v2; `pilot_lid_cad_results.md` v2;
> `headline_outcome_comparison.md` v3; `comparator_inspection.md`; `mvp_phe411_check.md`; `latent_ses_feasibility.md`
> v5; `plan_v0.6_review.md` v2.
>
> **Development history (drafts v0.2–v0.6, 2026-10-04 to 2026-10-07):** framing D5 option 1 (appraisal-led), then
> CAD-led with PAD, HF and stroke as pre-specified extensions of the appraisal; `qhet_mvmr` deferred (§10.3 B);
> mediation narrowed to an exploratory supplement (§9); the existing Xia 2025 gSES added as a narrow secondary analysis
> (§7.4).

**Working title:** *Can genetic instruments distinguish the effects of correlated socioeconomic measures on coronary
artery disease? An appraisal of multivariable Mendelian randomization designs, with extension to peripheral artery
disease, heart failure and stroke.*

---

## 0. Key design decisions (adopted)

| # | Decision | Adopted | Main reason | Changes an earlier decision? |
|---|---|---|---|---|
| D5 | Framing | Appraisal-led (option 1), **CAD-led**: can comparable SES-vs-SES models for CAD distinguish their components? PAD, HF and stroke are pre-specified extensions | CAD has three SES-vs-SES MVMR claims to appraise (Li D 2025; Le & Padmanabhan 2026 preprint, same EA4 and CAD files; Xu 2024 abstract); PAD has none. Expected SEs for PAD (MVP EUR) are ≈ 2.6× those for CAD (Aragam EUR) (`R/02i`). Chosen before any causal estimate | Yes: replaces the PAD-led emphasis of v0.2–v0.5 |
| C | Primary CAD file | **Aragam 2022 EUR (GCST90132314)** for all CAD analyses; Aragam multi-ancestry as sensitivity | Closer ancestry match: the exposure GWAS are predominantly European (primary TDI is the Pan-UKB multi-ancestry meta-analysis, so the EUR-only TDI sensitivity is retained), and the LD reference and LDSC scores are EUR; the Li D-comparable model uses Aragam EUR | Yes: v0.1–v0.5 had multi-ancestry as primary |
| D2/D3 | PAD outcome roles | **MVP `Phe_440_2` EUR is the primary PAD outcome** (PAD is an extension, not the headline). Analyses of the rounded van Zuydam data are **withheld from the main results** and reported only as **exploratory supplementary** analyses (§4.2). If full-precision values arrive before Phase 4, van Zuydam becomes the pre-specified secondary PAD outcome. Roles stay fixed whatever the results | MVP EUR is full-precision, PAD-specific, UK Biobank-free and has more cases (15,036 vs 12,086). The van Zuydam deposit is rounded (OR, SE and EAF to 2 decimals; p to 2 significant figures), and its SE conversion is provisional (§4.2) | **Yes: D2** had MVP as replication only |
| D7 | Sample overlap | Operational procedure in §2.1: a comparator without UKB participants for every outcome; LDSC intercepts from EUR-only files; Burgess-type calculations reported as illustrative scenarios, not as estimated bias. MRlap is not pre-specified | MRlap's assumption that overlap is independent of case/control status is not established. UK Biobank's share of controls is 2.3–9.2× its share of cases | Adds datasets (§17) |
| H | Harmonization | Adopt `HARMONIZATION_RULES_draft.md` draft v3 with the choices in §6 | Matching audit: 765–781 of 781 instruments are allele-consistent for each outcome | — |
| M | Mediation | **Secondary, exploratory supplement.** Education only, single-mediator, all four outcomes. Product of coefficients is the main estimand and the difference method secondary. CIs from a joint SNP bootstrap that keeps the dependence from shared instruments but not from shared participants. Strength gate as in §9.3. Mediator-specific proportions are not summed. MVP outcomes are not treated as overlap-reducing for the SBP and lipid models (§9.5) | Full-panel SES conditional F < 10; smoking fails in every single-mediator model; the independence assumptions of simple product/difference CIs do not hold | Narrows v0.1 H3 |
| R | "Robust finding" rule | **Removed.** The sensitivity battery characterizes uncertainty and is not a pass/fail gate (§8) | A binary label compounds multiple-testing and power problems | Replaces v0.1 §8 rule |
| E | MVMR estimators | **IVW-MVMR** (random effects) with its SEs, read alongside conditional F. **`qhet_mvmr` deferred** (not used in this paper). `mr_mvgmm` not used | `qhet_mvmr`: τ² search tolerance-limited (estimates moved by up to 0.075); df hard-coded (L − 2); its τ² profile searches negative values where its weights are invalid; final fits converged, but a strict every-step rule rejected all 30 test fits, leaving numerical questions unresolved (§10.3 B). Deferral is a scope decision | New |
| I | Appraisal interpretation | F < 10 (this study's reporting requirement) → weak-instrument-limited, with no attribution to one measure. F ≥ 10 does not by itself establish distinguishability or causal validity (§10.4) | Strength supports estimation; attribution also needs uncertainty and validity assumptions; methods designed for weak instruments can operate below 10 | Replaces the v0.3 rule |
| L | Latent SES factor | **Narrow secondary analysis** with the existing Xia 2025 gSES summary statistics (GCST90566708): univariable only, fixed "shared genetic SES liability" interpretation (§7.4). No factor GWAS is built here | Feasibility check (`latent_ses_feasibility.md` v5): open GRCh37 file; 379 instruments (mean F 50.6; relevance only); at most 1 instrument Q_SNP-heterogeneous; unconstrained rg with EA4 1.05 (SE 0.03) and with Kweon 0.99 (very high genetic similarity, not equivalence), so it cannot separate components | Yes: v0.5–v0.6 had deferred it; restores a narrowed v0.1 extension |
| F | FinnGen PAD endpoints | Keep `I9_ATHSCLE` as a broader-atherosclerosis sensitivity outcome. Drop `I9_OTHPER` and `I9_ARTEMBTHRLOW` | Under the DF12 definitions, none of the three is PAD-specific, and the last two are not atherosclerotic PAD | Narrows the PAD sensitivity set |

---

## 1. Research questions and hypotheses

### 1.1 Primary (appraisal, CAD)
- **Q1.** When correlated SES measures (educational attainment, income, occupational status, area deprivation) are
  modelled jointly in two-sample MVMR, how much conditional instrument strength is available to attribute CAD
  associations to an individual measure?
- **Q2.** Given that strength, how much uncertainty remains in IVW-MVMR estimates for individual SES measures across
  covariance assumptions and CAD data sources (Aragam EUR, Aragam multi-ancestry, MVP `Phe_411`)? Can comparable models
  of published SES-vs-SES CAD analyses distinguish their components?

### 1.2 Primary (application)
- **H1.** Higher genetically predicted SES (each of the four measures, univariable MR) is associated with lower risk of
  CAD, PAD, heart failure and any stroke (4 × 4 grid). The 16 cells form **one pre-specified primary testing family**.
  CAD leads the manuscript, and "extension" (PAD, HF, stroke) refers to the appraisal (§1.3, §10) only; it does not
  demote those outcomes within the H1 grid.

### 1.3 Secondary and extensions
- **Extension (PAD, HF, stroke).** The same appraisal (conditional F; IVW-MVMR with uncertainty) is reported for PAD,
  HF and any stroke under the same rules. Conditional F depends on the exposure instruments and is nearly identical
  across outcomes (instrument retention ≥ 99%); outcome precision differs.
- **Shared SES factor (secondary; §7.4).** Univariable MR of the Xia 2025 general SES factor (gSES) on CAD, read as the
  effect of genetically predicted shared SES liability. It does not identify which SES component acts and is not used
  to adjudicate between components.
- **H2 (exploratory).** The education–CVD association is consistent with partial mediation through established risk
  factors, assessed in single-mediator MVMR for education only (§9).
- **H3.** Univariable PAD associations are consistent across PAD definitions and data sources: MVP EUR, FinnGen broader
  atherosclerosis, and van Zuydam (exploratory unless full precision arrives). MVP African-ancestry PAD is analysed as
  exploratory transportability (§7.3).

### 1.4 Interpretation rules (pre-specified)
1. Findings from **reconstructed comparable models** (§10.2) describe those models only. They are never stated as
   properties of the published models they resemble, whose exact instrument sets and covariance assumptions are unknown.
   We do not claim that any published model had weak instruments.
2. The appraisal assesses what the instruments can distinguish and how much uncertainty remains. It does not claim that
   education is not causal, or that published estimates are wrong.
3. Conditional F = 10 is a **reporting requirement chosen for this study** and a diagnostic threshold. It is not a
   general condition for identification: methods designed for weak instruments can operate below 10, and strength
   alone does not establish causal validity. Each value is reported as a range across the
   four covariance scales (§10.1), with its distance from 10.
4. Write "highly genetically correlated", not "genetically identical".
5. "UK Biobank-free" means no shared UK Biobank participants. It does not mean overlap-free (§2.1).
6. Attenuation of an education estimate after adjusting for a mediator is not by itself evidence of mediation (§9.4).
7. Adequate conditional instrument strength supports estimation. It does not by itself show that the effects of
   correlated measures are distinguishable, causal, or attributable to one measure; those judgments also need the
   uncertainty of the estimates and the validity assumptions (§10.4).
8. **Scope of conclusions.** This paper can assess conditional instrument strength and illustrate conventional MVMR
   uncertainty. It cannot provide definitive corrected effects or determine which SES component truly dominates.

---

## 2. Study design

Two-sample univariable and multivariable MR using publicly available GWAS summary statistics; no individual-level data
are used. All exposure GWAS include UK Biobank (UKB). So does every primary outcome GWAS except PAD (MVP EUR). Overlap is
handled by design and quantified with approximate indicators (§2.1); no model-based correction is applied to the
primary estimates. A one-sample design is not used because no single cohort has adequately powered GWAS of all four SES
measures and all four outcomes.

### 2.1 Sample overlap (D7): operational procedure

**Step 1 — Comparator without UKB participants ("UKB-free") for every outcome.** These come from MVP gwPheWAS (dbGaP `phs002453`, open access, GIA
EUR stratum):

| Primary outcome | UKB-free comparator | Caveat |
|---|---|---|
| CAD (Aragam 2022 EUR) | MVP `Phe_411` ischaemic heart disease (EUR 141,465 cases / 275,809 controls; checked, `mvp_phe411_check.md`) | Broader EHR phecode in a predominantly male veteran population; expected SEs 1.11–1.22× Aragam EUR. Specification: rsID matching, SE from the 95% CI, no imputation exclusion (r² < 0.8 excluded as sensitivity) |
| PAD (MVP `Phe_440_2` EUR) | Primary is already UKB-free. The contrasts are van Zuydam (includes UKB) and FinnGen `I9_ATHSCLE` | — |
| Heart failure (HERMES 2.0) | MVP `Phe_428` heart failure (EUR 58,749 / 370,923; checked) | Phecode definition (nonhypertensive congestive HF); HERMES uses cohort-specific definitions. Expected SEs ≈ 1.5–1.6× HERMES |
| Any stroke (GIGASTROKE) | MVP `Phe_433_21` cerebral artery occlusion with infarction (EUR 20,766 cases) | **Approximates ischaemic stroke, not any stroke.** MVP has no any-stroke phecode, so this is not a like-for-like comparison (§4.3) |

MVP lists no cohort shared with any exposure GWAS. This is evidence class c: no shared cohort was found, which is not
proof that no individuals are shared. A UKB-free comparator therefore **reduces but does not eliminate** possible overlap,
and estimates from it also differ in phenotype definition and population (§4.4). FinnGen is not used as an overlap
comparator: it shares legacy cohorts with EA3, EA4 and Kweon (≤ 1.4% of exposure N) and contributes to HERMES 2.0 and
GIGASTROKE.

**Step 2 — Cross-trait LDSC intercepts** (exposure × outcome, 16 primary cells plus comparators). HapMap3 SNPs and
1000G EUR LD scores are used, so **outcome files must be EUR-only where one exists**:

| Outcome | File used for the LDSC intercept | Status |
|---|---|---|
| CAD | Aragam EUR (GCST90132314) | Downloaded, md5-verified |
| Any stroke | GIGASTROKE EUR (GCST90104539) | Downloaded 2026-10-06, md5-verified |
| Heart failure | HERMES 2.0 EUR-only | **Availability unconfirmed.** If no EUR-only file is public, the multi-ancestry file is used and the intercept is labelled *approximate* |
| PAD | MVP EUR; van Zuydam (EUR) | MVP: downloaded. van Zuydam: rounded values, so labelled *approximate* |
| Mediator models | GLGC lipids and DIAMANTE T2D are trans-ancestry | Intercepts labelled *approximate* |

The intercepts are reported as descriptive indicators of sample overlap and phenotypic correlation. They are not used
to correct univariable estimates.

**Step 3 — Burgess-type overlap scenarios (illustrative, not estimated bias).** Burgess, Davies & Thompson (2016,
doi:10.1002/gepi.21998) show that, under a null causal effect, the bias from sample overlap is approximately the
overlap fraction × 1/F, expressed relative to the observational association. We report this quantity for each
univariable cell as an **illustrative scenario**:

  scenario value = φ_UKB × (1/F),

where:
- **φ_UKB** is an upper-bound **UKB-mediated** overlap fraction: min(UKB N in the exposure GWAS, UKB N in the outcome
  GWAS) ÷ total outcome N, using the cohort counts documented in `d7_evidence.tsv`.
- **F** is the exposure's overall instrument F in its own selection GWAS (`INSTRUMENT_SUMMARY`): EA4 53.5, Kweon 45.8,
  TDI 34.3, occupation 41.5, Hill 39.8, ISEI 42.9. It is recomputed on each outcome-specific harmonized set before reporting.

The scenario values are **not** estimates of the bias in our results, for these reasons:
- the approximation is derived under a null causal effect;
- F is computed in the discovery sample of selected instruments and is likely overstated (winner's curse), which
  would understate the bias;
- φ_UKB counts only UK Biobank. It omits other cohorts shared between exposure and outcome GWAS, such as FinnGen
  legacy cohorts and other consortium members, and individuals shared without a shared cohort label;
- it ignores the case/control asymmetry of UKB's contribution (UKB's share of controls is 2.3–9.2× its share of
  cases), whose direction of effect on bias in a case-control outcome is not captured by the formula.

Where an exposure's UKB N is not documented (currently Kweon), the cell is "not calculated". For gSES (§7.4) no
scenario value is calculated: the factor has no single UKB N, and the UKB N of its largest component is not an
established upper bound on the factor's overlap bias. Its overlap is reported as unresolved. For MVP outcomes
**φ_UKB = 0**. That is the UKB-mediated overlap only, not total overlap, which is unknown (class c, Step 1).

Illustrative scenario values from current inputs, shown to document the calculation and not as reassurance about
bias. The first four columns use multi-ancestry outcome totals. The CAD (EUR) column uses UKB counts derived from the
documented shares (18.7% of 181,522 cases and 44.5% of 984,168 controls ≈ 471,900):

| φ_UKB × 1/F (% of observational association) | CAD (multi) | HF | Any stroke | PAD (van Zuydam) | CAD (EUR, headline) |
|---|---|---|---|---|---|
| EA4 (UKB 441,121) | 0.60 | 0.40 | 0.43 | 1.67 | 0.71 |
| TDI (UKB 439,198) | 0.93 | 0.62 | 0.67 | 2.61 | 1.10 |
| Occupation CAMSIS (UKB 273,157) | 0.48 | 0.33 | 0.41 | 1.43 | 0.56 |
| Hill income (UKB 286,301) | 0.52 | 0.37 | 0.45 | 1.56 | 0.62 |
| Kweon income | not calculated | not calculated | not calculated | not calculated | not calculated |

The main overlap evidence in the paper is the comparison with comparators that have no UKB participants (Step 1),
read with the phenotype and population differences in §4.3–4.4. The scenario values are supplementary.

**Step 4 — MRlap** is not pre-specified. Its overlap assumption is not established for any outcome.

**Step 5 — Overlap in MVMR.** The covariance-scaled conditional F (§10.1) assesses exposure-side strength only. The
estimators selected for this study do not correct exposure–outcome sample overlap. This is stated as a limitation, and MVMR results
on UKB-free comparators are reported alongside.

---

## 3. Exposures

| Exposure | GWAS | N | Role | Instruments (`R/02`) |
|---|---|---|---|---|
| Educational attainment | Okbay 2022 EA4, 23andMe excluded | 765,283 | primary | 408 |
| Income | Kweon 2025 Income Factor | 668,288 (effective) | primary | 131 |
| Area deprivation (TDI) | Pan-UKB 189-irnt (meta) | 439,198 | primary | 9 |
| Occupational status | Akimova 2024 CAMSIS (GCST90446160) | 273,157 | primary | 63 |
| Education (comparator) | Lee 2018 EA3 | 766,345 | sensitivity | 306 |
| Income (comparator) | Hill 2019 household income (UKB only) | 286,301 | sensitivity; Li D-comparable model | 27 (Li D model: S1 set, §10.2) |
| TDI (comparator) | Pan-UKB 189-irnt EUR-only stratum | ≈ 421,945 (estimated) | sensitivity (rg with meta = 1.04) | 8 |
| Occupation (comparator) | Akimova 2024 ISEI (GCST90446163) | 271,769 | Le-comparable model only | 52 |
| Shared SES factor | Xia 2025 gSES, Genomic SEM common factor of EA, household income, CAMSIS and TDI (GCST90566708; UKB + SSGAC) | 947,466 (effective) | secondary (§7.4) | 379 |

Units are SD per allele. TDI is sign-flipped so that higher means less deprived. All exposures include UKB. gSES effects
are reported per released factor unit (higher = higher SES). The units of the released SNP effects are not stated in the main
text, the supplementary information or the deposit metadata. Supplementary Table 4 is consistent with a unit-variance genetic-factor
parameterization. However, the relationship between that parameterization and the released SNP-effect units remains
unresolved. Comparison with our EA4 associations does not establish a specific rescaling (§7.4, §18).

---

## 4. Outcomes

### 4.1 Primary outcomes and comparators

| Outcome | Primary GWAS | Cases | Sensitivity / comparators |
|---|---|---|---|
| **CAD (headline)** | **Aragam 2022 EUR (GCST90132314)** | 181,522 | Aragam multi-ancestry (GCST90132315; 210,842 cases); MVP `Phe_411` ischaemic heart disease (UKB-free; 141,465 EUR cases) |
| **PAD** | **MVP `Phe_440_2` EUR** | 15,036 | FinnGen R12 `I9_ATHSCLE` (broader atherosclerosis: ICD-10 I70, ICD-9 440/4370/4371X, ICD-8 440); van Zuydam 2021 (§4.2); MVP AFR (exploratory) |
| Heart failure | HERMES 2.0 all-cause (GCST90726617) | 153,174 | MVP `Phe_428`; HERMES 2.0 EUR-only if available |
| Any stroke | GIGASTROKE multi-ancestry (GCST90104534) | 110,182 | GIGASTROKE EUR (GCST90104539); MVP `Phe_433_21` |

Outcome subtypes (HERMES 2.0 HF subtypes, GIGASTROKE ischaemic subtypes) are **not in the committed scope**: their files
were not downloaded or screened before the freeze. The phenotype-mismatch caveats for the comparators (§4.3) stand. A
subtype analysis could be added only by a dated amendment logged in `results/logs/DEVIATIONS.md` before any of its results
are examined; one added after examining results would be reported as post hoc.

### 4.2 van Zuydam 2021 PAD (D3 rule)
- The public deposit (figshare 12723404) gives OR, `or_se` and EAF rounded to two decimals and p rounded to two
  significant figures. No full-precision or UKB-excluded copy is public (checked: figshare, GWAS Catalog, the CVD
  Knowledge Portal, and the supplement). The authors were contacted on 2026-10-04.
- **If full-precision values arrive before Phase 4:** van Zuydam is the pre-specified secondary PAD outcome, analysed
  like the other outcomes.
- **Otherwise:** van Zuydam analyses are **withheld from the main results** and reported only in the supplement as
  *exploratory, rounded-data* analyses, under `HARMONIZATION_RULES_draft.md` draft v3, rule 11:
  - β = ln(OR as deposited).
  - **SE(β) = `or_se` / OR — supported but provisional.** Evidence (`R/02g`; variant-level table
    `results/tables/vanzuydam_se_scale_check.tsv`): over the joint rounding range of OR and p, in a fixed set of 16
    large-effect variants, 12 are compatible only with `or_se` being SE(OR), 1 only with SE(log OR), and 3 with
    neither. In a wider pre-stated set of 1,045 variants the counts are 346 / 15, with 521 compatible with both and
    163 with neither. The deposited p-values are systematically less extreme than OR and `or_se` imply: the median
    ratio of `or_se` to the p-implied SE(OR) is 0.94–0.95. The cause is unknown.
  - **Sensitivity:** SE = max(`or_se`/OR, |ln OR| / z(p)), the more conservative of the two.
  - **Variants with OR = 1.00 are retained**, with β = 0 and their SE. No variant is dropped for rounding, and no sign
    is assigned or imputed. Currently 16–33% of matched instruments have OR = 1.00 (rsID matching before harmonization).
  - The proportion of instruments with OR = 1.00, and the rounding intervals, are reported with every estimate.
- Roles do not change in response to results.

### 4.3 Stroke definition caveat
The primary outcome is **any stroke** (GIGASTROKE AS). The UKB-free comparator `Phe_433_21` approximates **ischaemic
stroke**. Differences between the two can reflect phenotype (haemorrhagic strokes are included only in the primary
outcome), as well as overlap or population. They are not interpreted as overlap bias.

### 4.4 MVP limitations
MVP gwPheWAS outcomes are EHR-derived phecodes in a predominantly male US veteran population with a distinctive
socioeconomic and exposure profile (military service, smoking prevalence, VA care access). Phecode case definitions
differ from consortium definitions, and controls exclude related phecodes. Differences between MVP and consortium
estimates can arise from phenotype, population and selection, not only from overlap. MVP results are reported on
the log-OR scale, with SE derived from the reported 95% CI (§6).

---

## 5. Mediators

| Mediator | GWAS | N | Ancestry | Units |
|---|---|---|---|---|
| BMI | Yengo 2018 (GIANT + UKB) | ~700,000 | EUR | SD |
| Lifetime smoking index | Wootton 2020 | 462,690 | EUR (UKB) | SD |
| LDL-C, HDL-C, triglycerides | Graham 2021 (GLGC) | 1.65 M | trans-ancestry | SD (inverse-normal) |
| Systolic blood pressure | Keaton 2024 (UKB + ICBP + MVP + BioVU) | ~1.03 M | EUR | mmHg |
| Type 2 diabetes | Mahajan 2022 (DIAMANTE) | 180,834 cases / 1,159,055 | trans-ancestry | log-OR |

Use is restricted as described in §9. Smoking and the trans-ancestry mediators carry the caveats in §9.3 and §2.1.
Keaton SBP and GLGC lipids include MVP data. All seven mediator GWAS include UKB; for DIAMANTE 2022 (Supplementary Tables 1–2, checked 2026-10-07), UKB
contributes 19,119 cases / 423,698 controls (EUR) and 1,404 / 6,849 (SAS). MVP is not among DIAMANTE's 84 listed
studies: no shared cohort was found (class c, §2.1), which does not exclude shared individuals (§9.5).

---

## 6. Instrument selection and harmonization

**Selection** (`R/02_select_instruments.R`, run): p < 5×10⁻⁸; clumping r² < 0.001 within 10 Mb (1000G EUR); F > 10;
ambiguous palindromes with MAF > 0.42 removed.

**Harmonization** (`HARMONIZATION_RULES_draft.md` draft v3):

| Rule | Choice |
|---|---|
| Matching | GRCh37 outcomes: chr:pos + allele pair, with an rsID cross-check. GRCh38 (MVP): rsID + allele pair. No liftover |
| Palindromes | TwoSampleMR `action = 2` (frequency-based orientation). MVP AFR: all palindromes dropped |
| Missing instruments | No proxies in primary analyses (≤ 2% lost per outcome). Proxies with r² ≥ 0.8 as a sensitivity analysis |
| MVP imputation quality | No exclusion (no instrument has r² < 0.3). Sensitivity: exclude r² < 0.8 |
| MVP effect scale | β = ln(OR); SE = (ln upper − ln lower) / (2 × 1.959964), validated against reported p. For the GWAS Catalog copies (`Phe_428`, `Phe_433_21`) the `ci_lower`/`ci_upper` columns are used; their `standard_error` column is not |
| van Zuydam | §4.2 rules |
| MR-Steiger | Directionality reported for every cell. Filtering is a sensitivity analysis only |
| MHC | No exclusion. Sensitivity: exclude chr6:25–34 Mb |
| Residual LD | `plink --r2` audit reported. Sensitivity: drop the weaker SNP of each pair with r² ≥ 0.001 |
| Retention | Reported per cell and per step in `results/tables/harmonization_retention.tsv` |

MVMR instrument sets are the union of the component exposures' instruments, re-clumped at r² < 0.001 and harmonized
per outcome. Conditional F is **recomputed on each outcome-specific harmonized set**; the current values (§9, §10) are
exposure-side estimates on the pre-harmonization union.

---

## 7. Univariable analysis (application, H1)

### 7.1 Primary
The estimator is IVW with multiplicative random effects. The grid is 4 exposures × 4 outcomes = 16 cells, tested at
Bonferroni α = 0.05/16 = 0.003125. Estimates are reported as OR per SD of exposure, with 95% CIs.

### 7.2 Comparators
Lee EA3, Hill income, EUR-only TDI, Aragam multi-ancestry, GIGASTROKE EUR, the MVP comparators, FinnGen `I9_ATHSCLE`, and van Zuydam (exploratory
unless full precision arrives) are reported as sensitivity analyses without multiplicity adjustment.

### 7.3 Exploratory transportability
MVP AFR `Phe_440_2` with EUR-derived instruments, all palindromes dropped. This is exploratory and is not a replication:
the instruments, LD and allele frequencies were derived in EUR samples.

### 7.4 Shared SES factor (secondary)
- **Exposure:** Xia 2025 gSES (GCST90566708), instruments selected by the rules of §6 (379 instruments; mean F 50.6;
  R² 0.020). No instrument is excluded on Q_SNP heterogeneity (at most one would qualify).
- **Outcomes:** CAD Aragam EUR (main), with MVP `Phe_411` EUR as the comparator without UK Biobank participants and Aragam multi-ancestry as
  sensitivity. PAD (MVP EUR), HF (HERMES 2.0) and any stroke (GIGASTROKE) follow only as part of the pre-specified
  extensions, with their §17 comparators.
- **Estimators:** IVW (multiplicative random effects) with the §8 sensitivity battery. Mean F 50.6 supports instrument
  relevance only; pleiotropy, overlap and the interpretation of a composite exposure are assessed separately.
- **Overlap:** UKB contributes to every gSES indicator, so overlap with UKB-containing outcomes is qualitatively similar to
  EA4's but not demonstrably identical (composition differs; the effective N blends indicators). EA4's overlap fractions
  are not transferred, and no numerical §2.1 scenario is calculated for gSES: its overlap is reported as unresolved and
  discussed qualitatively. MVP `Phe_411` is a comparator without UK
  Biobank participants, not an independent check: unknown participant overlap remains possible, and both analyses use
  the same gSES exposure estimates.
- **Interpretation (fixed):** "genetically predicted shared SES liability", reported as a potentially precise estimate
  under the stated MR assumptions. The estimate is not attributed to any component and is not an intervention effect.
  Given very high rg with EA4 and Kweon, similarity to their estimates is likely but not guaranteed; neither agreement
  nor disagreement is evidence about which component acts.
- **Units:** per released gSES factor unit. The units of the released SNP effects are not stated in the main text, the
  supplementary information or the deposit metadata. Supplementary Table 4 is consistent with a unit-variance genetic-factor
  parameterization. However, the relationship between that parameterization and the released SNP-effect units remains
  unresolved. Comparison with our EA4 associations does not establish a specific rescaling. We therefore report effects
  per released factor unit.
  Released units are preserved, and no effect magnitude is compared across exposures.
- **Not done:** MVMR of gSES with any of its components (near-collinear by construction); comparison of
  magnitudes with component estimates (unit unresolved); the UKB-only factor (GCST90566707).

---

## 8. Sensitivity battery (every univariable estimate)

| Method | What it characterizes |
|---|---|
| MR-Egger (slope, intercept) | Directional pleiotropy under InSIDE |
| Weighted median | Estimate valid if < 50% of weight comes from invalid instruments |
| Weighted mode | Estimate valid under the plurality-valid assumption |
| MR-PRESSO | Global heterogeneity; outlier-corrected estimate |
| Cochran's Q, I² | Heterogeneity |
| Leave-one-out; single-SNP | Influence of individual instruments |
| MR-Steiger | Directionality under its assumptions; not confirmation |
| MR-RAPS | Many weak instruments and balanced pleiotropy; not an overlap correction |

**Use (replaces the v0.1 "robust finding" rule).** The battery characterizes uncertainty and is **not a pass/fail
gate**. No estimate is labelled "robust" or "not robust". For each cell, the report states:
- the direction and magnitude of each estimator, with its CI;
- the Egger intercept and its CI;
- heterogeneity (Q, I²) and the PRESSO global test, with the outlier-corrected estimate;
- the largest change in the leave-one-out estimate, and which SNP produced it.

Agreement between estimators with different assumptions is described as such. Divergence is described along with the
assumption that would explain it. Primary inference rests on the IVW estimate and its pre-specified α (§12).

---

## 9. Mediation (H2): education only, exploratory

### 9.1 Scope
- **Status: secondary, exploratory supplement.** Mediation results are reported in the supplement with the
  limitations in §9.4–9.5. They do not enter the paper's primary conclusions or abstract claims.
- **Exposure:** EA4 only. No mediation estimates are produced for income, TDI or occupation. Their conditional F is
  < 10 in the full mediator panel and in most single-mediator models.
- **Mediators:** each of the seven in §5, in a **single-mediator** model (EA + one mediator).
- **Outcomes:** all four primary outcomes (CAD = Aragam EUR; PAD = MVP EUR). Mediation is run regardless of the
  univariable results. It is not conditional on any robustness gate.

### 9.2 Estimands
- **Main — two-step product of coefficients.** Indirect effect = β(EA→M) × β(M→Y | EA), where:
  - β(EA→M) is the univariable IVW estimate of EA on the mediator, using EA instruments;
  - β(M→Y | EA) is the mediator coefficient from the single-mediator IVW-MVMR.
- **Secondary — difference method.** Total (univariable EA→Y) minus direct (EA coefficient in the MVMR).
- The proportion mediated is reported only when the total effect's 95% CI excludes 0, and then with its CI.
- **Proportions for different mediators are not summed or combined.** Single-mediator models are fitted separately,
  and their pathways can overlap (e.g. BMI → T2D → outcome), so a sum has no interpretation as a joint proportion.

### 9.3 Strength gate
A single-mediator model's estimates are **eligible for interpretation** only if **both** EA and the mediator have
conditional F ≥ 10 at **all four covariance scales** (§10.1), in the outcome-specific harmonized instrument set.
Passing the gate does not by itself make an estimate interpretable as mediation (§9.4). Models failing the gate are
reported in the supplement, labelled *weak-instrument-limited*, and not interpreted.

Current exposure-side values (`R/02d`, pre-harmonization union) will be **reassessed after harmonization**:
- EA conditional F is 11.5–35.9 across single-mediator models.
- Mediator conditional F: BMI 36.8–42.4, HDL 129–136, LDL 150–151, SBP 61–62, T2D 45, TG 126–130.
- **Smoking is 6.55–7.30**, so the education–smoking models are currently expected to be weak-instrument-limited.

The full seven-mediator panel is reported only as a strength diagnostic. At the LDSC-intercept scale the SES
conditional F is EA 6.6, income 2.9, occupation 2.5 and TDI 1.3 (1.2–7.1 across all scales); smoking is 2.2–3.2.

### 9.4 Uncertainty
The two steps share instruments, and the GWAS share participants, so the estimates are not independent. Simple
formulas that assume independence (product method) or ignore the covariance between total and direct effects
(difference method) are therefore not used as the main intervals.
- **Main CI — joint SNP bootstrap.** Instruments for the model (the union of EA and mediator instruments, harmonized
  for the outcome) are resampled with replacement as one set, 2,000 times. In each resample all component estimates are
  recomputed from that resample: β(EA→M), the MVMR coefficients, and the total effect, each using the instruments of
  its own type that fall in the resample. Product and difference estimates are then formed. Percentile 95% intervals
  are reported. This keeps the dependence that comes from shared and overlapping instruments.
- **Failed resamples.** A resample counts as failed if any component estimate errors or is non-finite, or if the
  resample contains too few instruments of a type to fit it (fewer than 3). With any failure, the interval is labelled
  **failure-affected** and exploratory, because failed resamples can differ systematically from successful ones.
  Above 5% failures the interval is **withheld**.
- **What it does not capture:** correlation in estimation errors due to shared participants between the EA, mediator
  and outcome GWAS (EA4 and all seven mediator GWAS include UKB; DIAMANTE 2022 includes UKB, §5). The
  estimators selected for this study do not correct mediation estimates for this. It is stated as a limitation, and the intervals are described as not accounting for
  sample overlap. Joint resampling of instruments captures dependence from shared instruments only; the intervals are
  **not** presented as fully accounting for dependence between the estimates, and the mediation analyses remain
  exploratory for this reason as well.
- For comparison, the delta-method interval assuming independence is shown in the supplement, labelled as such.

### 9.5 Overlap in mediation models
- **MVP is not overlap-reducing for every mediator.** The SBP GWAS (Keaton 2024) includes MVP. The lipid GWAS
  (GLGC, Graham 2021) includes MVP data (MVP is a consortium author). DIAMANTE 2022 T2D lists no MVP cohort among its 84
  studies (checked 2026-10-07): no shared cohort found (class c), so MVP overlap is possible but not documented. BMI
  (Yengo 2018) and smoking (Wootton 2020) contain no MVP data.
- Mediation models with an MVP outcome and an SBP or lipid mediator are therefore labelled *mediator–outcome overlap:
  MVP*, and are not presented as overlap-reduced comparisons. Models with an MVP outcome and the T2D mediator are
  labelled *mediator–outcome overlap: MVP possible (no shared cohort found)*, and are likewise not presented as
  overlap-reduced.

### 9.6 Interpretation caveats
- **Noncollapsibility.** Binary outcomes are analysed on the log-OR scale, where total ≠ direct + indirect in general.
  Product and difference estimates can therefore diverge even without bias, and neither is an exact decomposition. For
  the low-prevalence outcomes the discrepancy is expected to be small, but it is not assumed to be zero.
- **Attenuation is not demonstrated mediation.** A smaller EA coefficient after mediator adjustment can also reflect
  pleiotropy, mediator measurement error, weak conditional instruments, or collider effects. Results are described as
  "consistent with mediation through X", not as precise proportions.
- SBP is in mmHg, T2D in log-OR, and the other mediators in SD. Indirect effects are reported on the outcome log-OR
  scale per SD of EA.

---

## 10. Appraisal component (Q1, Q2)

### 10.1 Conditional instrument strength
Sanderson–Windmeijer conditional F is computed for every MVMR model. Exposure-error covariance is specified at **four
scales: 0, 0.5, 1 and 1.5 × the LDSC-intercept-based correlation ρ** (pairwise ρ from cross-trait LDSC intercepts
between exposures). Results are reported as a range across the four scales. Implementation: `R/02b`, independently
verified by `R/02c` (96/96 values including the Le-comparable models, max relative difference 1.3×10⁻¹⁴). Recomputed on each outcome-specific harmonized
set.

### 10.2 Models

**A. Project models** (CAD first: Aragam EUR, Aragam multi-ancestry and MVP `Phe_411`; then the PAD, HF and stroke
extensions with their comparators):
1. Four-SES joint model: EA4 + Kweon income + TDI + Akimova CAMSIS.
2. EA4 + Kweon income.
3. EA4 + TDI.
4. EA4 + CAMSIS.
5. Hill variant: model 1 with Hill income in place of Kweon.

**B. Comparable models** (reconstructions; `reconstruction_feasibility.md` v2). Exact reproduction is not established
for any of them, and the interpretation rule in §1.4(1) applies.
1. **Li D 2025 CAD:** EA4 + Hill income, instrument set S1 (407 SNPs; `R/02e`), outcome Aragam EUR. Exposure-side
   conditional F: EA 2.24–3.06, income 2.03–2.55. LD-pruned sensitivity (383 SNPs): EA 2.20–2.99, income 1.97–2.46. The
   diagnostic set S2 (433 SNPs) fails the LD audit and is not estimated.
2. **Le 2026 single-adjustment models:** EA4 + Hill income (407 SNPs; conditional F EA 2.24–3.06, income 2.03–2.55);
   EA4 + Akimova ISEI (403 SNPs; EA 1.11–1.73, ISEI 1.06–1.58; ISEI downloaded 2026-10-06, 52 instruments)
   (`R/02b`; arithmetic verified on saved inputs by `R/02c`). These are additional comparable models, not
   verification of Le's models. EA4 + Hill uses the same 407-SNP set as the Li D-comparable model, so the two
   reconstructions share one instrument-strength result. Le's joint
   four-measure model is not reconstructed because its composition is uncertain.
3. **Ye 2025 HF (optional):** EA3 / income / TDI as OpenGWAS IDs `ieu-a-1239`, `ukb-b-7408`, `ukb-b-10011`, with HERMES 2020
   (`ebi-a-GCST009541`). Included only if the files are accessible and approved; otherwise reported as not reconstructed.

Published SNP counts are compared as checks only; no rule is tuned to match them.

### 10.3 Estimators

**A. Core: IVW-MVMR.** `MendelianRandomization::mr_mvivw` (multiplicative random effects), for every model in §10.2.
Each coefficient is reported with its SE and 95% CI, next to its conditional F range (§10.1), and on both the
primary outcomes and the UKB-free comparators. With weak conditional instruments, IVW-MVMR estimates are biased
and ordinary confidence intervals may be unreliable. **The direction of the bias is model-dependent** when exposures
are correlated: it is not necessarily toward zero, and it can differ between coefficients in the same model (Zhu,
Burgess & Grant, arXiv 2203.08668). In our `R/02f` simulation (one design, conditional F 2.4) the means were 0.210 and
−0.012 against true values of 0.3 and −0.1, i.e. attenuation in that design. That is an illustration, not a general
direction. Coefficients with conditional F < 10 are therefore labelled
weak-instrument-limited (§10.4) and not interpreted as effect sizes. No other MVMR estimator enters the main results;
`mr_mvgmm` is not used.

**B. `MVMR::qhet_mvmr` is deferred** (version 0.4.8; Q-statistic minimisation intended to be robust to weak
instruments). It is not used in this paper. The reasons, from the package source and tests on simulated data
(`R/utils/qhet_checks.R`, `R/02h_qhet_checks_test.R`, `results/tables/qhet_checks_test_*.tsv`), are:
- **τ² is tolerance-limited.** The package profiles a heterogeneity term τ² with `optimize` on (−10, 10) at the
  default tolerance (≈ 1.2×10⁻⁴, the same order as SNP–outcome variances). In 30/30 simulated datasets it returned
  τ² = ±2.28×10⁻⁵. Re-solving τ² at tolerance 10⁻¹⁰ moved the estimates by a median (maximum) of 0.023 (0.033) at
  conditional F ≈ 2.4, 0.005 (0.006) at F ≈ 12–21, and 0.027 (0.075) in a three-exposure model.
- **The τ² search enters invalid regions.** At negative trial values the weight denominators collapse to about 0, and
  the inner Nelder–Mead fits stop with a degenerate simplex (code 10). This happened in 4–7 of 28 profile evaluations
  in every dataset. Those trial points are discarded by the search, and the **final fits converged in 30/30 datasets**.
  A strict rule that every trial evaluation must succeed rejected all 30 fits (and every bootstrap resample). That rule
  is stronger than requiring a valid final solution, so **0/30 passing it does not show that the 30 final estimates were
  invalid**. It shows that the implementation's search behaviour would need further work before its results could be
  checked to the standard set here.

- **The df target is hard-coded** as L − 2 for any number of exposures.
- **Why defer rather than use with a weaker rule.** The tolerance-limited τ², the hard-coded df, the invalid search
  region and unestablished CI coverage together leave unresolved numerical and inference questions. Resolving them would
  turn this paper into a package-repair project. Choosing a weaker rule after seeing the strict one fail would also be
  results-driven. Deferral is therefore a scope decision, not a finding that the estimator is wrong.
- **What was verified is implementation only.** An independent re-implementation reproduced the package's estimates
  exactly, the final-fit objective agreed with the best of 12 starts (largest relative gap 1.2×10⁻⁷), a forced invalid
  τ² was correctly flagged, and bootstrap failure counting worked. None of this establishes statistical
  correctness or CI coverage.
- In a best-case simulation (`R/02f`: no pleiotropy or overlap), the package estimates were close to the truth at
  F ≈ 2.4. That suggests the approach may be useful, but repairing and validating the implementation is outside the
  scope of this paper. A weak-instrument-robust comparison is listed as future work (§16).

### 10.4 Interpretation rules (appraisal)
1. Conditional F is a diagnostic of strength, reported as a range across the four covariance scales (§10.1).
2. **F ≥ 10 at all four scales is the reporting requirement chosen for this study** before an IVW-MVMR coefficient
   is interpreted as more than weak-instrument-limited. It is a study convention, not a general condition for
   identification.
3. A coefficient with conditional F < 10 at any scale is labelled **weak-instrument-limited**. No statement ranks SES
   measures, or attributes the association to one measure, on the basis of such coefficients. For these models the
   appraisal reports the strength and the width of the uncertainty, not an attribution.
4. Meeting the requirement does not by itself establish that the effects are distinguishable, causal or attributable to
   one measure. Such a statement must also rest on the coefficient's CI, consistency across covariance scales and
   between the primary outcomes and UKB-free comparators, and the validity assumptions (pleiotropy diagnostics;
   overlap, §2.1). Untested assumptions are named.
5. Statements about comparable models follow §1.4(1).

---

## 11. Power and minimum detectable effect

**Method (computed before any causal estimate; `R/02i_outcome_precision.R`).** For each univariable cell, the expected
SE of the IVW estimate is 1/√Σ(β_X²/SE_Y²) over the instruments available in the outcome. It uses exposure effects from
the discovery GWAS and outcome SEs only. The minimum detectable OR per SD at 80% power is exp(±(z₁₋α/2 + z₀.₈) × SE),
with α = 0.003125 for the primary grid and α = 0.05 for comparators. Values are optimistic (winner's curse in β_X) and are
recomputed on each outcome-specific harmonized set. For gSES (§7.4) the same benchmark is computed from the released
betas and is expressed per released factor unit; it is not a standardized effect.

**Current values (α = 0.003125; OR per SD giving 80% power, protective / harmful):**

| Exposure | CAD Aragam EUR | PAD MVP EUR | HF HERMES 2.0 | Any stroke GIGASTROKE (multi-ancestry, primary) |
|---|---|---|---|---|
| EA4 | 0.93 / 1.08 | 0.82 / 1.21 | 0.93 / 1.07 | 0.92 / 1.09 |
| Kweon income | 0.87 / 1.16 | 0.69 / 1.45 | 0.87 / 1.15 | 0.85 / 1.18 |
| CAMSIS occupation | 0.87 / 1.15 | 0.70 / 1.43 | 0.88 / 1.14 | 0.85 / 1.17 |
| TDI | 0.59 / 1.69 | 0.25 / 4.00 | 0.62 / 1.61 | 0.56 / 1.79 |

**Pre-specified power statements.**
- These are the effect sizes giving 80% power under the stated assumptions. Smaller effects can still be detected, with
  lower probability. No cell is declared uninformative in advance; every estimate is reported and interpreted with its
  interval.
- **PAD:** income and occupation have 80% power near OR 0.69–0.70 per SD. Whether that is adequate depends on the
  relevant effect size (for reference, the published education → CHD estimate was 0.67 per SD; Tillmann 2017).
  **Deprivation (TDI) has limited power to detect modest effects in every outcome, most of all PAD** (80% power only
  near OR 0.25 / 4.0 per SD).
- A non-significant estimate in a cell with limited power is reported as compatible with the effects inside its
  confidence interval, not as evidence of no effect.
- **MVMR:** the conventional SE of each IVW-MVMR coefficient, √diag((BᵀWB)⁻¹), is reported as a precision benchmark
  across outcomes (four-SES model, CAD EUR: 0.067–0.114 per coefficient; PAD MVP EUR: 0.173–0.293). It treats exposure
  associations as known, so it reflects neither weak-instrument bias nor overlap or pleiotropy, and it is not a bound
  on uncertainty.

## 12. Multiple testing
- Primary univariable grid: one pre-specified family of 16 cells (4 exposures × CAD, PAD, HF, any stroke), Bonferroni
  α = 0.003125. CAD-led presentation does not change this family.
- Appraisal (§10): estimation-focused. Results are reported with CIs and conditional F; no hypothesis test drives a
  conclusion.
- Shared SES factor (§7.4), mediation, comparators and exploratory analyses: unadjusted, interpreted in the context of the primary correction.

## 13. Reporting and registration
- STROBE-MR checklist in the supplement.
- Code: R scripts in the GitHub repository, archived on Zenodo with a DOI at submission. `renv.lock` pins packages, and
  each run log records sessionInfo.
- Data: summary statistics referenced by DOI/accession and checksum in `data/raw/MANIFEST.tsv`. Nothing is redistributed
  unless the licence permits.
- Registration: this plan, version-frozen, is registered on OSF before any Phase 4 estimate is computed.

## 14. Deviations
Any deviation during execution, such as an inaccessible dataset or a method failing for a cell, is logged in
`results/logs/DEVIATIONS.md` with the deviation, reason, date and downstream effect. The log is published with the
paper. Adopting this plan logs each change from v0.1 there.

## 15. Manuscript outline
1. Introduction: education-primacy claims in SES→CAD MR (Li D 2025; Le & Padmanabhan 2026; Xu 2024); correlated
   SES measures.
2. Methods.
3. Results:
   1. Genetic correlations and the conditional strength of SES MVMR.
   2. Comparable-model reconstructions and IVW-MVMR uncertainty.
   3. Univariable SES→CVD: the 4 × 4 grid, with CAD first; the shared SES factor (§7.4) as a secondary estimate.
   4. Extensions: the appraisal repeated for PAD (MVP EUR), HF and stroke.
   5. Education mediation (single-mediator; supplement).
4. Discussion: what genetic instruments can and cannot distinguish, and implications for interpreting published claims.

## 16. Out of scope
Within-family MR; bidirectional MR; building or modifying a latent SES factor GWAS, and MVMR of a factor with its
components (the existing gSES is used univariably only, §7.4); PheWAS; ancestry-matched exposure analyses (no
adequately powered non-European SES GWAS); model-based overlap correction of MVMR estimates; weak-instrument-robust
MVMR estimation (`qhet_mvmr` deferred, §10.3 B; future work).

## 17. Additional datasets: status

| Dataset | Use | Condition |
|---|---|---|
| MVP `Phe_411` (EUR) | CAD UKB-free comparator | **Downloaded and checked 2026-10-06 (user-approved; `mvp_phe411_check.md`):** 141,465 cases / 275,809 controls; 778/781 instruments allele-consistent; expected SEs 1.11–1.22× Aragam EUR. Ischaemic heart disease is broader than Aragam CAD, so agreement or disagreement cannot be attributed solely to overlap |
| MVP `Phe_428` (EUR) | HF UKB-free comparator | **Downloaded and checked 2026-10-07** (GWAS Catalog GCST90475981, identical to the dbGaP release; `comparator_inspection.md`): 58,749 / 370,923; 778/781 instruments; SEs ≈ 1.5–1.6× HERMES 2.0 |
| MVP `Phe_433_21` (EUR) | Stroke UKB-free comparator (ischaemic approximation) | **Downloaded and checked 2026-10-07** (GWAS Catalog GCST90477991): 20,766 / 420,219; 778/781 instruments; SEs ≈ 1.6–1.7× GIGASTROKE EUR |
| GIGASTROKE EUR (GCST90104539) | EUR LDSC intercepts; stroke sensitivity | **Downloaded and checked 2026-10-06** (`comparator_inspection.md`): 649/657 instruments; SEs ≈ 1.19× multi |
| Akimova ISEI (GCST90446163) | Le-comparable CAD model (EA4 + ISEI) | **Downloaded and checked 2026-10-06:** 52 instruments; positions GRCh37 despite the metadata label |
| Xia 2025 gSES (GCST90566708) | Shared SES factor, secondary (§7.4) | **Downloaded and checked 2026-10-07** (`latent_ses_feasibility.md` v5): md5 matches EBI; 379 instruments; LDSC h² 0.101 (paper 9.4%); unconstrained rg EA4 1.05, Kweon 0.99, CAMSIS 0.93, TDI 0.58. Paper's accession list truncated, so the release is identified from Catalog metadata; released SNP-effect units unresolved (not stated in the paper, supplement or deposit metadata) |
| HERMES 2.0 EUR-only | EUR LDSC intercepts | Only if publicly available |
| OpenGWAS `ieu-a-1239`, `ukb-b-7408`, `ukb-b-10011`; HERMES 2020 `ebi-a-GCST009541` | Ye-comparable HF model (**optional**; extension) | Only if accessible; omission is not a deviation |

## 18. Open items at the freeze (none blocks registration)
- The van Zuydam reply (full precision; possibly a UKB-free GoLEAD version). It determines the §4.2 branch.
- The Kweon authors' reply (UKB-excluded Income Factor), and Kweon's UKB N for §2.1 Step 3.
- Accessibility of the Ye OpenGWAS files.
- Optional author query (Xia/Hill): "What normalization and units were used for the released gSES SNP effects, and how
  do they relate to the factor parameterization in Supplementary Table 4? Were any transformations applied before
  release?" Not blocking: §7.4 reports per released factor unit and compares no magnitudes.
- Accepted as limitations, not to be resolved before registration: IVW-MVMR weak-instrument bias (handled by labelling,
  §10.4); the absence of an overlap correction for MVMR and mediation estimates; no weak-instrument-robust MVMR
  estimator (`qhet_mvmr` deferred, §10.3 B).
- None of the open author replies (van Zuydam, Kweon) blocks registration. The van Zuydam branch in §4.2 is
  pre-specified either way.
- **Other outcomes considered and not added:** venous thromboembolism and abdominal aortic aneurysm (distinct vascular
  comparisons; novelty and precision unchecked), atrial fibrillation (crowded; Kaisaier 2025). Outcome subtypes
  (GIGASTROKE ischaemic, HERMES 2.0 HF) are outside the committed scope (§4.1). Adding an outcome would not address the difficulty of
  separating correlated SES exposures, which is the paper's subject.
