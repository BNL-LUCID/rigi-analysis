# File Format Reference

Column schemas at each pipeline step. Pickle and CSV outputs share the
same columns at any given step.

## 1. `*_seqinfo.txt` (SigProfiler output, step 1)

Tab-delimited, one row per mutation per sample. Format produced by
`SigProfilerMatrixGenerator`.

| Column   | Example          | Notes |
|----------|------------------|-------|
| Sample   | `d0_W0`          | `<dose>_<timepoint>` |
| Type     | `SNV`            | Mutation class |
| Position | `12345`          | 1-based genomic coordinate |
| Context  | `N[C>T]`         | Trinucleotide context with substitution |
| Ref      | `C`              | Reference allele |
| Alt      | `T`              | Alternate allele |
| Strand   | `0`              | Strand orientation |

Files are emitted per chromosome under
`<output>/<TYPE>/.../<chr>_seqinfo.txt`.

## 2. `{chr}_{type}_mutations.pkl` / `.csv` (preprocessing output, step 2)

| Column                     | Description |
|----------------------------|-------------|
| `Chromosome`, `Start`, `End` | Genomic coordinates (Start is 0-based) |
| `Sample`                   | `<dose>_<timepoint>` |
| `Ref`, `Alt`               | Reference and alternate alleles |
| `Context`                  | Original SigProfiler context string |
| `Mutation_Strand`          | Strand orientation |
| `Dose`                     | `Control`, `dA`–`dE` |
| `Timepoint`                | `W0`–`W3` |
| `Mutation_Type`            | `SNV`, `DBS`, `MNS`, `ID` |
| `MutationID`               | Per-row identifier: `chr_pos_ref_alt_sample` |
| `Quality_Annotation`       | SigProfiler quality flag (`N` / `Q`) |
| `Transcription_Annotation` | Transcribed strand flag (`T` / `U`) for SNV/DBS/MNS |
| `Indel_Type`               | `Ins` / `Del` (ID only) |
| `Indel_Mechanism`          | `M` / `C` / `R` (ID only) |
| `Indel_Size`               | Indel length in bp (ID only) |

## 3. `all_<TYPE>_annotated.pkl` / `.csv` (annotation output, step 4)

All columns from step 2, plus:

| Column            | Description |
|-------------------|-------------|
| `Gene_Name`       | Overlapping gene symbol (`Unknown` for intergenic) |
| `Feature_Type`    | `exon`, `intron`, `5utr`, `3utr`, `promoter`, `Intergenic` |
| `Gene_Strand`     | `+` / `-` / `NA` |
| `Gene_Location`   | `Genic` / `Intergenic` |

The fixed output filename produced by `mutation_annotation.py` is
`annotated_mutations.pkl`. Step 5 expects `all_<TYPE>_annotated.pkl`,
so rename after each annotation run.

## 4. `mutation_annotations_dose_<dose>_<TYPE>.csv` (categorical patterns, step 5)

One row per unique mutation (deduplicated across timepoints).

| Column                   | Description |
|--------------------------|-------------|
| `MutationID`             | `chr_pos_ref_alt` (no sample suffix; cross-timepoint key) |
| `Chromosome`, `Start`, `Ref`, `Alt` | Genomic coordinates and alleles |
| `W0`, `W1`, `W2`, `W3`   | Categorical state per timepoint (e.g., `W1_Treatment`, `W2_Lost`) |
| `Pattern`                | Compact 4-character string (e.g., `0T00`) |
| `Category`               | Biological grouping (e.g., `Treatment_Only_W1`) |

Companion `all_mutations_dose_<dose>_<TYPE>.csv` carries 0/1 binary
presence across all sample-timepoint combinations.

## 5. `<TYPE>_dose_<dose>_merged.csv` (merged output, step 6)

Inner join of the annotated pickle (step 4) and the categorical pattern
CSV (step 5) on `chr_pos_ref_alt`. Contains every column from both,
plus:

| Column          | Description |
|-----------------|-------------|
| `Pattern_Group` | High-level grouping: `Radiation-specific`, `Control-specific`, `Baseline`, `Other` |

This is the file consumed by the SV-side scripts for breakpoint-proximal
enrichment and dose-stratified analyses.

## 6. Sankey JSONs (step 7)

### `all_chromosomes_trajectories.json` (preferred)

Keys are full 4-week trajectory strings; values are mutation counts.

```json
{
  "W0_Absent->W1_Exposed->W2_Lost->W3_Exposed_Recurrent": 4744,
  "W0_Absent->W1_Lost->W2_Lost->W3_Lost":                 89221,
  ...
}
```

Each mutation contributes to exactly one trajectory.

### `all_chromosomes_flows.json` (legacy / sanity check)

Pairwise transitions, computed by collapsing trajectories. Two-state
keys, useful only for cross-checking trajectory totals.

```json
{
  "W2_Lost->W3_Exposed_Recurrent": 5759,
  "W2_Lost->W3_Lost":              95281,
  ...
}
```
