# Conventions (SV side)

Conventions used throughout the SV pipeline. The mutation pipeline at
`../final codes_mutation/` uses overlapping conventions; pattern string
notation (`PATTERNS.md`) is shared.

## Sample naming

- **`d0`** — unirradiated control
- **`dA`–`dE`** — radiation-exposed samples at experimentally determined
  dose rates of **0.36, 0.20, 0.40, 1.47, 2.62 mGy/hr** respectively
- **Timepoints** — `W1`, `W2`, `W3` (weekly post-exposure); `W0` is the
  pre-exposure baseline
- **Filename convention** — `d0_vs_d<X>_W<n>_annotated.tsv` after AnnotSV

## Pattern strings

Four-character trajectory `[W0][W1][W2][W3]`:

- **`T`** — present in treated samples
- **`C`** — present in control samples
- **`B`** — present in both
- **`0`** — absent

W0 is always treated as `0` in the leading position (pre-exposure
baseline; SV calling pre-exposure is not meaningful for radiation
attribution).

The seven radiation-specific SV patterns:
`0T00`, `00T0`, `000T`, `0TT0`, `0T0T`, `00TT`, `0TTT`.

See [`../mutation/patterns.md`](../mutation/patterns.md) for the full pattern
encoding (shared with the mutation side).

## Dose bins

For dose-stratified analyses (Fig 5, §3.5):

- **Low** — `dA`, `dB`, `dC` (0.20–0.40 mGy/hr)
- **High** — `dD`, `dE` (1.47–2.62 mGy/hr)

The threshold (1 mGy/hr) separates dose rates well below typical
occupational exposure limits from those substantially above them.

## AnnotSV `Annotation_mode == 'full'`

All SV-counting scripts filter to `Annotation_mode == 'full'` before
counting. AnnotSV emits one row per overlapping gene by default (so a
SV touching N genes appears N times); the `'full'` rows are the
per-SV summary rows that should be counted exactly once.

If you forget this filter, every count downstream (SV totals, repeat
involvement percentages, breakpoint enrichments) is inflated by a
factor proportional to mean genes-per-SV.

## Inversion deduplication

AnnotSV occasionally reports the same INV with reversed coordinates
(start/end swapped) due to how it consolidates breakend (BND) calls
back into INV events. Scripts that count INVs collapse these by using
sorted `(start, end)` coordinate keys.

## Breakpoint matching

- **Cross-timepoint SV tracking** — 1000 bp tolerance on breakpoint
  positions, with `SV_type` and `Chromosome` exact match required.
  Accommodates positional jitter from independent variant calling.
- **Mutation-to-breakpoint windows** — 10, 25, 50, 100 bp. The 10 bp
  window is the headline analysis (manuscript §3.4); larger windows
  document distance decay.
- **Mutation tracking across timepoints** — exact coordinate match
  (`chr_pos_ref_alt`). No tolerance applied to point mutations.
