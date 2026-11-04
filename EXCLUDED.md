# Files omitted from the public release

The working repository is private. This public release was assembled from it by `release/make_public_snapshot.py`
(kept in the private repository). Files were omitted for one of three reasons. All omitted data files are regenerated
by the pipeline from the source downloads.

## 1. Confirmed restrictions or conditions

| Omitted | Reason |
|---|---|
| `data/instruments/ea_okbay2022_instruments.tsv`, `ea_lee2018_instruments.tsv` | SSGAC terms prohibit republication of the summary statistics. SSGAC permitted (2026-10-11) the effect allele, beta and SE of the instrument SNPs only; these files also carry allele frequency, P, N, F and R². The permitted columns are in `manuscript/tables/instrument_lists.tsv`. |
| `data/instruments/income_kweon2025_instruments.tsv` and Kweon-derived per-variant inputs | The Kweon et al income-factor summary statistics are licensed CC BY-NC-SA; redistribution is permitted only under those conditions. Omitted; per-variant files derived from them would have to carry the same licence. |

## 2. Source terms not verified (omitted conservatively)

Per-variant exposure, mediator and outcome effect sizes from the other sources (Pan-UKB, Akimova et al, Xia et al,
Hill et al, GIANT, Wootton et al, GLGC, Keaton et al, DIAMANTE, CARDIoGRAMplusC4D/Aragam et al, MVP, HERMES, GIGASTROKE,
FinnGen, van Zuydam et al). Their redistribution terms were not checked for this release; omission does not imply that
redistribution is prohibited.

- `data/instruments/*_instruments.tsv` (all exposure and mediator instrument files)
- `results/tables/condF_inputs/`, `results/tables/condF_inputs_mediators/`, `results/tables/pilot_lid_cad_inputs/`
  (per-variant beta and SE matrices used for conditional F)
- `results/tables/vanzuydam_se_scale_check.tsv` / `.xlsx` (per-variant van Zuydam outcome odds ratios)
- `results/reference/*.rds` (cell-level reference objects containing harmonized per-variant data)
- `data/cell_inputs/outcome_n_per_snp.tsv.gz` (per-variant outcome sample sizes)

## 3. Derived per-variant outputs from which source effects could be reconstructed

Combined with published exposure effects (for example the EA instrument betas in eTable 19), these would allow outcome
effects to be recovered.

- `results/tables/mr_univariable_singlesnp.tsv` (per-SNP ratio estimates and SEs)
- `results/tables/mr_univariable_loo.tsv` (leave-one-out estimates for every SNP)
- `results/tables/mr_univariable_steiger_snp.tsv` (per-SNP exposure and outcome variance explained)

## Also not included

Working notes, correspondence, internal reviews, manuscript drafts, journal-specific submission files, literature
files and other working files. The manuscript itself is published by the journal.
