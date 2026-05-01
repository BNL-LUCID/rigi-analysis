# Complete Mutation Analysis Pipeline: VCF to Final Analysis

## Pipeline Overview

```
VCF Files
    ↓
[1] sigprofiler.py → *_seqinfo.txt files (per chromosome, per mutation type)
    ↓
[2] mutation_preprocessing.py → {chr}_{type}_mutations.pkl
    ↓
[3] annotation_preprocessing.py → hg38_interval_trees.pkl + hg38_feature_data.pkl
    ↓
[4] mutation_annotation.py → all_{type}_annotated.pkl
    (parallel_annotate.py wraps this for SNV; see Parallel helpers section)
    ↓
[5] mutation_pattern_assignment.py →
    ├─ mutation_annotations_dose_{dose}_{type}.csv (C/T/B/0 patterns)
    └─ all_mutations_dose_{dose}_{type}.csv (binary 0/1 patterns)
    (parallel_pattern_assignment.py runs doses concurrently; see Parallel helpers)
    ↓
[6] merge_annotation.py → {type}_dose_{dose}_merged.csv
    ↓
[7] compute_sankey.py  → dose_{dose}/all_chromosomes_trajectories.json (per-dose, QC)
                         + combined/all_chromosomes_trajectories.json (manuscript view)
                         + matching all_chromosomes_flows.json files (legacy pairwise format)
    render_sankey.py   → publication-quality PNG from a trajectories JSON
    ↓
[8] Analysis scripts → Use merged files for downstream analysis
```

---

## Environment Setup

The pipeline runs on Linux + Python 3.9. Step 1 (`sigprofiler.py`) optionally
uses a CUDA GPU via PyTorch; all other steps are CPU-only.


```bash
# 1. New conda env, Python 3.9
conda create -n lucid_mut python=3.9 -y
conda activate lucid_mut

# 2. Install all Python deps (numpy, pandas, scipy, matplotlib, seaborn,
#    intervaltree, requests, SigProfilerMatrixGenerator, SigProfilerExtractor)
pip install -r requirements.txt

```

### CPU-only machine

The default `pip install -r requirements.txt` pulls a
CPU-only PyTorch as a transitive dep of SigProfilerExtractor; that's enough.
Run `sigprofiler.py` without `--gpu`.

### Reference genome (offline workflow)

`sigprofiler.py` normally downloads `GRCh38.tar.gz` from the AlexandrovLab FTP
on first run. On VMs without outbound FTP (SciServer, some HPC nodes), pre-stage
the tarball locally and pass `--offline-genome /path/to/GRCh38.tar.gz`.

The script copies the tarball into the SigProfilerMatrixGenerator install
tree, extracts it, and **auto-flattens nested layouts** — so it handles both
the flat `<genome>/{1.txt,2.txt,...,MT.txt}` layout and the wrapped
`<genome>/chromosomes/tsb/<genome>/...` layout that `tar -czf` from a parent
directory produces. On the second invocation it sees the genome is already
extracted and short-circuits.

---

## Detailed Step-by-Step Guide

### Step 1: SigProfiler Matrix Generation
**Script:** `sigprofiler.py`

**Input:**
- VCF files (directory containing `.vcf` or `.vcf.gz` files)
- Optional: pre-downloaded `GRCh38.tar.gz` for offline genome install

**Output:**
- `{chromosome}_seqinfo.txt` files for each mutation type (SNV, DBS, ID)
- Format: Tab-delimited text with mutation information

**Command (full — GPU + offline genome):**
```bash
python sigprofiler.py \
    --input            ./vcf_files \
    --output           ./sigprofiler_output \
    --reference        GRCh38 \
    --project          HUVEC_radiation \
    --min-signatures   1 \
    --max-signatures   10 \
    --nmf-replicates   100 \
    --cpu              -1 \ ##use all available CPU, makes the first matrix generation step fast##
    --gpu \
    --offline-genome   /path/to/GRCh38.tar.gz
```

**Common variations:**
```bash
# CPU-only (no GPU available)
python sigprofiler.py -i ./vcf_files -o ./sigprofiler_output --offline-genome /path/to/GRCh38.tar.gz

# Pin a specific GPU on a multi-GPU box
CUDA_VISIBLE_DEVICES=0 python sigprofiler.py -i ./vcf_files -o ./sigprofiler_output --gpu --offline-genome /path/to/GRCh38.tar.gz

# Run in background, log to file
nohup python sigprofiler.py -i ./vcf_files -o ./sigprofiler_output --gpu --offline-genome /path/to/GRCh38.tar.gz > sig.log 2>&1 &
tail -f sig.log

```

**Argument reference:**

| Flag | Default | Purpose |
|------|---------|---------|
| `-i / --input` | required | Directory of `.vcf` / `.vcf.gz` files |
| `-o / --output` | required | Output directory for SigProfiler results |
| `-r / --reference` | `GRCh38` | Reference genome (`GRCh37` or `GRCh38`) |
| `-p / --project` | `SigProfilerProject` | Project name (used in output filenames) |
| `-m / --min-signatures` | `1` | Lower bound of NMF signature search |
| `-M / --max-signatures` | `10` | Upper bound of NMF signature search |
| `-n / --nmf-replicates` | `100` | NMF replicates per signature count |
| `-c / --cpu` | `-1` | CPU count (`-1` = all available) |
| `--gpu` | off | Run NMF on CUDA GPU via PyTorch |
| `--offline-genome` | unset | Local `{reference}.tar.gz` to bypass FTP |



**What it does:**
- Processes VCF files through SigProfilerExtractor (which invokes
  SigProfilerMatrixGenerator internally to build the per-chromosome
  trinucleotide matrices and emit `*_seqinfo.txt` files)
- Generates mutation matrices by type (SNV, DBS, ID, MNS)
- Performs de novo NMF signature extraction (CPU or GPU)
On a 72 core machine with P100 GPU, the code takes about 2.5 hours for 15 vcf files. 
---

### Step 2: Mutation Preprocessing
**Script:** `mutation_preprocessing.py`

**Input:**
- `*_seqinfo.txt` files from Step 1
- Expected to live under `<input-dir>/{SNV,DBS,ID,MNS}/` — searched
  recursively, so SigProfiler's nested layout works directly

**Output:**
- `{chromosome}_{type}_mutations.pkl` - Structured mutation data
- `{chromosome}_{type}_mutations.csv` - Human-readable version
- `{chromosome}_{type}_summary.csv` - Quality control summary

**Where SigProfiler writes its output (important):**

`SigProfilerMatrixGenerator` always writes back into the input VCF directory,
producing a tree like:

```
vcf_files/                       ← original VCF dir (--input to sigprofiler.py)
└── output/
    └── vcf_files/               ← SigProfiler uses the basename as project name
        ├── SNV/
        │   └── ... <chr>_seqinfo.txt   (possibly nested another level deep)
        ├── DBS/
        ├── ID/
        └── MNS/
```

`mutation_preprocessing.py` searches each `<TYPE>/` subdirectory **recursively**,
so it doesn't matter how deep SigProfiler nests the files — point `--input-dir`
at the directory that contains the four type subfolders.

**Command:**
```bash
# Point --input-dir at SigProfiler's output tree
python mutation_preprocessing.py \
    --input-dir   ./vcf_files/output/vcf_files \
    --chromosome  all \
    --output      ./processed_data \
    --summary     ./summary_data

# Process a single chromosome (useful for debugging)
python mutation_preprocessing.py \
    --input-dir   ./vcf_files/output/vcf_files \
    --chromosome  chr1 \
    --output      ./processed_data \
    --summary     ./summary_data

# Default: --input-dir is '.', so if you cd into the SigProfiler output
# directory first you can omit the flag entirely
cd ./vcf_files/output/vcf_files
python /path/to/mutation_preprocessing.py -c all -o ./processed_data -s ./summary_data
```

**Argument reference:**

| Flag | Default | Purpose |
|------|---------|---------|
| `-i / --input-dir` | `.` | Directory containing `SNV/`, `DBS/`, `ID/`, `MNS/` (typically `<vcf_dir>/output/<vcf_dir_basename>/`) |
| `-c / --chromosome` | `all` | Specific chromosome or `all` |
| `-o / --output` | `processed_data` | Output directory for `*_mutations.pkl` / `.csv` |
| `-s / --summary` | `summary_data` | Output directory for QC summary CSVs |
| `-t / --mutation-types` | `SNV DBS MNS ID` | Subset of types to process. Useful when SigProfiler has finished some types but not others, or to run two instances in parallel. |

**Splitting the run** (recommended when SigProfiler hasn't finished all types):

`SigProfilerMatrixGenerator` runs ID as a separate pass after SNV/DBS/MNS, so
those three are typically ready well before ID. You don't need to wait — kick
off preprocessing for the ready types now and finish ID later:

```bash
python mutation_preprocessing.py \
    --input-dir       ./vcf_files/output/vcf_files \
    --output          ./processed_data \
    --summary         ./summary_data \
    --mutation-types  SNV DBS MNS

# Later, after SigProfiler finishes ID:
python mutation_preprocessing.py \
    --input-dir       ./vcf_files/output/vcf_files \
    --output          ./processed_data \
    --summary         ./summary_data \
    --mutation-types  ID
```

Outputs accumulate into the same `processed_data/` directory across runs,
so no merging is needed.

**Running types in parallel.** The script processes types sequentially within
one process. If you want true parallelism, run two instances pointing at the
same `--input-dir` and `--output` but different `--mutation-types`:

```bash
python mutation_preprocessing.py -i ... -o processed_data -t SNV DBS &
python mutation_preprocessing.py -i ... -o processed_data -t MNS ID &
wait
```

(Modest gain — parsing is mostly I/O-bound — but worthwhile on fast disks.)

**What it does:**
- Recursively finds and parses `*_seqinfo.txt` under each requested type subfolder
- Extracts: chromosome, position, ref, alt, sample, dose, timepoint
- Adds mutation context (quality, transcription, indel mechanism)
- Creates permanent mutation IDs: `chr_pos_ref_alt_sample`

**Key columns created:**
- `Chromosome`, `Start`, `End`
- `Sample`, `Ref`, `Alt`
- `Dose`, `Timepoint`
- `Mutation_Type` (SNV/DBS/ID)
- `MutationID`
- `Quality_Annotation`, `Transcription_Annotation`
- `Indel_Type`, `Indel_Mechanism`, `Indel_Size`

---

### Step 3: Annotation Preprocessing (One-time setup)
**Script:** `annotation_preprocessing.py`

**Input:**
- UCSC refGene database (auto-downloaded)

**Output:**
- `hg38_interval_trees.pkl` - Genomic interval trees
- `hg38_feature_data.pkl` - Gene/feature annotations

**Command:**
```bash
python annotation_preprocessing.py \
    --build hg38 \
    --annotation-dir annotations
```

**What it does:**
- Downloads UCSC refGene annotations
- Builds interval trees for fast genomic lookups
- Identifies: genes, exons, introns, UTRs, promoters
- One-time setup, reuse for all analyses

---

### Step 4: Mutation Annotation
**Script:** `mutation_annotation.py`
**Parallel wrapper:** `parallel_annotate.py` — recommended for SNV (see
[Parallel helpers](#parallel-helpers))


**Input:**
- A mutations pickle from Step 2. Either an `all_<TYPE>_mutations.pkl`
  (combined-chromosomes; recommended) or per-chromosome pickles.
- `hg38_interval_trees.pkl` and `hg38_feature_data.pkl` from Step 3

**Output (per `-o` directory):**
- `annotated_mutations.pkl` / `.csv` — annotated mutations (single file)
- `quality_transcription_analysis/` — quality and transcription QC plots/CSVs
- `indel_mechanism_analysis/` — indel size/type/frameshift breakdowns (only
  produced for `Mutation_Type == ID`; SNV/DBS/MNS skip this section cleanly)

**Command (combined-chromosomes input, recommended):**
```bash
# Run once per mutation type. Each type writes to its own -o directory so
# fixed output filename `annotated_mutations.pkl` doesn't get overwritten.
for mut_type in SNV DBS ID MNS; do
    python mutation_annotation.py \
        -m processed_data/all_${mut_type}_mutations.pkl \
        -a annotations \
        -b hg38 \
        -o annotated_mutations/${mut_type}
done

# Step 5 expects filenames matching `*_<TYPE>_annotated.pkl` in a flat dir,
# so rename the fixed output for each type after Step 4 completes:
for mut_type in SNV DBS ID MNS; do
    src="annotated_mutations/${mut_type}/annotated_mutations.pkl"
    [ -f "$src" ] && mv "$src" "annotated_mutations/${mut_type}/all_${mut_type}_annotated.pkl"
    src_csv="annotated_mutations/${mut_type}/annotated_mutations.csv"
    [ -f "$src_csv" ] && mv "$src_csv" "annotated_mutations/${mut_type}/all_${mut_type}_annotated.csv"
done
```

**Command (per-chromosome input):**
```bash
python mutation_annotation.py \
    -m processed_data/chr1_SNV_mutations.pkl \
    -a annotations \
    -b hg38 \
    -o annotated_mutations/chr1_SNV
```

**What it does:**
- Annotates each mutation with genomic features
- Identifies gene names, feature types (exon/intron/UTR)
- Calculates distances to nearest genes
- Adds strand information

**Key columns added:**
- `Gene_Name`
- `Feature_Type` (exon, intron, 5utr, 3utr, promoter)
- `Gene_Strand`
- `Gene_Location`

---

### Step 5: Pattern Assignment
**Script:** `mutation_pattern_assignment.py`
**Parallel wrapper:** `parallel_pattern_assignment.py` — runs all doses
concurrently (see [Parallel helpers](#parallel-helpers))


**Input:**
- `*_<TYPE>_annotated.pkl` files from Step 4 (e.g., `all_DBS_annotated.pkl`).
  The script globs the input dir for files matching that suffix and
  concatenates them. With the combined-file convention from Step 4 there's
  one file per type.

**Output (per dose found in the data):**
- **Categorical patterns:** `<pattern_dir>/<TYPE>/dose_<dose>/mutation_annotations_dose_<dose>_<TYPE>.csv` (and `.pkl`)
  - Columns: `MutationID, Chromosome, Start, Ref, Alt, W0, W1, W2, W3, Pattern, Category`
  - `MutationID` is `chr_pos_ref_alt` (the deduplicated, cross-timepoint key)
- **Binary patterns:** `<pattern_dir>/<TYPE>/dose_<dose>/all_mutations_dose_<dose>_<TYPE>.csv` (and `.pkl`)
  - 0/1 presence across all sample-timepoint combinations

**Command:**
```bash
# Run once per mutation type — emits one categorical + one binary file per
# dose found in the data (e.g. dose_dA/, dose_dB/, ...).
for mut_type in SNV DBS ID MNS; do
    python mutation_pattern_assignment.py \
        -i annotated_mutations/${mut_type} \
        -o pattern_analysis/${mut_type} \
        -m ${mut_type}
done
```

**Note on `--by-chromosome`:** the flag derives "chromosomes" from filename
prefixes (the part before the first `_`). With a single combined input
`all_<TYPE>_annotated.pkl` it sees one prefix (`"all"`), loops once, and
behaves identically to the default mode — safe but redundant. The flag only
becomes hazardous if the input dir contains **both** a combined file and
per-chromosome files: each iteration writes to the same fixed output
filename per dose, so later iterations overwrite earlier ones. With purely
combined-file input, you can ignore the flag.

**What it does:**
- Groups mutations by PermanentMutationID (chr_pos_ref_alt)
- Determines presence in control vs treatment at each timepoint
- Assigns categorical states:
  - `W0_Present` or `W0_Absent` (baseline)
  - `W1/W2/W3_Control`, `W1/W2/W3_Treatment`, `W1/W2/W3_Both`, `W1/W2/W3_Lost`
- Creates compact pattern strings (C/T/B/0)
- Categorizes into biological groups

**Output format (categorical):**
```csv
MutationID,W0,W1,W2,W3,Category,Pattern,Chromosome,Start,Ref,Alt
chr1_100_A_T,W0_Absent,W1_Treatment,W2_Lost,W3_Lost,Treatment_Only_W1,0T00,chr1,100,A,T
```

---

### Step 6: Merge Annotations with Patterns
**Script:** `merge_annotation.py`

**Input:**
- Annotated file from Step 4 (combined-chromosomes preferred):
  `<annotated-dir>/<TYPE>/all_<TYPE>_annotated.csv`
- Pattern file from Step 5:
  `<pattern-dir>/<TYPE>/dose_<dose>/mutation_annotations_dose_<dose>_<TYPE>.csv`

The script auto-detects combined-file vs. per-chromosome layouts; combined
mode runs as a single join, per-chromosome mode loops chromosomes.

**Output:**
- `<output-dir>/merged_data/<TYPE>_dose_<dose>_merged.csv` — merged dataset
  (annotations + Pattern + Pattern_Group + W0–W3 states)
- Per-dose plots under `<output-dir>/<TYPE>/dose_<dose>/`: pattern × gene
  location, pattern × feature type, strand bias, top genes, etc.

**Command:**
```bash
# Dose values must match the directory names emitted by Step 5 (e.g. dA,
# not A, since pattern_analysis/<TYPE>/dose_dA/...).
python merge_annotation.py \
    --annotated-dir annotated_mutations \
    --pattern-dir pattern_analysis \
    --output-dir merged_analysis \
    --mutation-types SNV DBS ID MNS \
    --doses dA dB dC dD dE
```

**What it does:**
- Joins annotations to patterns on `chr_pos_ref_alt` (inner join — annotated
  mutations not present in any temporal pattern are dropped, expected behaviour)
- Adds `Pattern_Group`: Radiation-specific / Control-specific / Baseline (CBBB) / Other
- Generates summary plots per (type, dose)
- With `--compare-doses`, also produces cross-dose comparison plots in
  `<output-dir>/<TYPE>/dose_comparison/`

**Output contains:**
- All genomic features (Gene_Name, Feature_Type, Gene_Location, Gene_Strand)
- Temporal patterns (W0, W1, W2, W3 states + Pattern string + Category)
- Pattern_Group label (Radiation-specific / Control-specific / Baseline / Other)
- Sample, Dose, Timepoint context

---

### Step 7: Sankey Visualization
**Scripts:** `compute_sankey.py` (trajectories) + `render_sankey.py` (figure)

Two-step process: first compute trajectories from annotated mutations, then
render the figure. Operates directly on Step 4 annotated pickles (not Step 5
patterns) because the manuscript state taxonomy needs the full
per-(sample, timepoint) presence matrix.

**Why trajectories, not pairwise transitions:** the renderer is
*trajectory-aware* — it draws each mutation's full 4-week path as a
contiguous ribbon. An earlier draft used per-transition pairwise flows
(e.g. `"W1_Exposed->W2_Lost": 1234`), but stacking ribbons independently at
each node created visual artifacts where two unrelated ribbons could appear
to flow continuously through a path that no mutation actually followed. The
trajectory format (`"W0_Absent->W1_Exposed->W2_Lost->W3_Exposed_Recurrent": 1234`)
preserves end-to-end identity, so the visual matches the underlying data.

**State taxonomy** (matches manuscript Fig 2):
- W0: `Present`, `Absent`
- W1: `Both`, `Lost`, `Control`, `Exposed` (no Recurrent — W0 is baseline)
- W2: `Both`, `Exposed_Recurrent`, `Control_Recurrent`, `Exposed`, `Control`, `Lost`
- W3: same as W2

`*_Recurrent` fires only when the mutation was present in that arm at some
prior timepoint **and** absent in that arm at the immediately preceding
week. Mutations newly appearing without prior presence are classified as
`Exposed`/`Control`, not Recurrent.

**Inputs/Outputs:**
- Input: `annotated_mutations/<TYPE>/all_<TYPE>_annotated.pkl` (Step 4)
- Output (per `--output-dir`):
  - `dose_<dose>/all_chromosomes_trajectories.json` — per-dose, QC
  - `combined/all_chromosomes_trajectories.json` — single-pass full dataset (manuscript view)
  - Matching `all_chromosomes_flows.json` (pairwise transitions, kept for legacy/sanity-check use)

**Why combined and not per-dose for the manuscript figure:** at this state
taxonomy the temporal dynamics are largely dose-agnostic (per-dose Sankeys
look essentially identical), so the manuscript shows one combined figure.
Dose-specific signal is captured separately in fig 5 / dose-stratified analyses.

**Why not sum per-dose JSONs to get combined:** every mutation present only
in controls appears in *all 5* per-dose subsets, so summing would
over-count control-side flows by 5×. `compute_sankey.py` produces the
combined JSON directly from the full dataset to avoid this.

**Command:**
```bash
# Step 7a: compute trajectories (writes both per-dose and combined JSONs)
python compute_sankey.py \
    --input      annotated_mutations/DBS/all_DBS_annotated.pkl \
    --output-dir sankey_flows/DBS

# Step 7b: render the combined manuscript figure
python render_sankey.py \
    --trajectories-json sankey_flows/DBS/combined/all_chromosomes_trajectories.json \
    --output            sankey_figures/DBS_combined.png \
    --title             "Temporal Dynamics for DBS" \
    --subtitle          "All Doses, All Chromosomes" \
    --dpi               300

# Optional: render an individual dose for QC
python render_sankey.py \
    --trajectories-json sankey_flows/DBS/dose_dA/all_chromosomes_trajectories.json \
    --output            sankey_figures/DBS_dose_dA.png \
    --title             "Temporal Dynamics for DBS" \
    --subtitle          "Dose dA, All Chromosomes"
```

**Sanity check after computing:**
```bash
# Total trajectory counts == total unique mutations (counted once)
python -c "
import json
t = json.load(open('sankey_flows/DBS/combined/all_chromosomes_trajectories.json'))
print(f'unique trajectories: {len(t):,}')
print(f'total mutations:     {sum(t.values()):,}')
"
```

**Diagram features:**
- Muted-tone palette matching the manuscript (Both / Exposed / Control / Lost / Recurrent / Baseline)
- Trajectory-aware bezier ribbons: each mutation's full 4-week path drawn as one contiguous unit
- Width auto-scales to mutation count; small trajectories visible alongside dominant ones
- 16:9 publication-ready PNG at requested DPI

---

### Step 8: Downstream Analysis
**Scripts:** Various analysis scripts

**Input:**
- `{type}_dose_{dose}_merged.csv` from Step 6

**Common analyses:**
- Differential gene analysis
- Pathway enrichment
- Hotspot identification
- Temporal pattern analysis
- Dose-response curves
- Mutation signature analysis

---

## Parallel helpers

Two outer-loop wrappers sit alongside the core scripts and parallelize the
two slowest steps without modifying the underlying logic. They emit the same
output layout as the sequential scripts — drop-in replacements when the
sequential runs are too slow.

### `parallel_annotate.py` (Step 4 wrapper)

Chunk-level parallelism for `mutation_annotation.py`. Splits the input
`all_<TYPE>_mutations.pkl` into N row-balanced chunks, runs N annotation
subprocesses concurrently, and concatenates the per-chunk annotated outputs
into a single combined `all_<TYPE>_annotated.pkl`.

**When to use:** SNV (10s of M rows). Brings ~6h sequential down to ~30-40min
on a multi-core box. DBS/MNS/ID are usually fast enough sequentially.

**Command:**
```bash
python parallel_annotate.py \
    --input            processed_data/all_SNV_mutations.pkl \
    --annotation-dir   annotations \
    --build            hg38 \
    --output           annotated_mutations/SNV/all_SNV_annotated.pkl \
    --workers          24 \
    --no-csv           # optional: skip the (large) CSV write for huge tables
```

**Notes:**
- `--output` must be a `.pkl` file path, not a directory. The CSV is written
  alongside (filename derived by replacing `.pkl` → `.csv`) unless `--no-csv`.
- Each worker loads its own copy of the interval-tree annotation data
  (a few hundred MB to ~1 GB on hg38). RAM scales linearly with `--workers` —48
- One failed chunk fails the whole run with exit 1. Per-chunk logs land in
  `<output>_parallel_work/out/chunk_NNN/subprocess.log`. The scratch dir is
  preserved on failure (so you can recover via `pd.concat` of the surviving
  chunk pickles); deleted on success unless `--keep-work-dir` is passed.

### `parallel_pattern_assignment.py` (Step 5 wrapper)

Dose-level parallelism for `mutation_pattern_assignment.py`. Splits the
annotated input by dose (each split contains the controls + that one dose),
launches one subprocess per non-control dose, each writing to its own
`dose_<dose>/` subdirectory. No merge step — each subprocess produces its
own final per-dose output.

**When to use:** any type where sequential pattern assignment is too slow.
Capped at the number of non-control doses in the data (typically 5).

**Command:**
```bash
python parallel_pattern_assignment.py \
    --input          annotated_mutations/SNV/all_SNV_annotated.pkl \
    --output-dir     pattern_analysis/SNV \
    --mutation-type  SNV \
    --workers        5
```

**Notes:**
- Output layout is identical to a sequential run. You can verify on a small
  type first by diffing against an existing sequential output:
  ```bash
  diff <(sort old/dose_dA/mutation_annotations_dose_dA_DBS.csv) \
       <(sort new/dose_dA/mutation_annotations_dose_dA_DBS.csv)
  # expected: no output
  ```
- Memory caveat: each worker loads its filtered DataFrame (controls + one
  dose ≈ ~30% of full input). 5 workers × ~30% = ~1.5× peak vs. sequential.
  On large inputs (SNV) this can OOM-kill workers — kernel sends SIGKILL,
  workers die with returncode -9 and empty logs (no chance to flush). If you
  see this pattern, drop `--workers` to 2 or 3, or run failed doses
  sequentially using the per-dose split inputs in the work dir (preserved
  on failure).
- The script captures each subprocess's output into
  `<output-dir>_parallel_work/logs/<dose>.log`. Empty logs after a SIGKILL
  is the OOM-kill signature.

### Caveats common to both

- The helpers don't modify the underlying scripts, so all the bug fixes and
  guards in `mutation_annotation.py` / `mutation_pattern_assignment.py` apply
  inside each subprocess.
- Per-chunk QC plots from `mutation_annotation.py` (the
  `quality_transcription_analysis/` and `indel_mechanism_analysis/`
  directories) are produced per-chunk and discarded with the scratch dir —
  they're not representative when computed on a slice. If you want global QC
  plots after parallel annotation, run `mutation_annotation.py` once on the
  merged output (skipping reannotation by passing the already-annotated file
  as `-m` won't work — the QC code is inside `main()` after annotation).
  In practice the merged DataFrame is what matters; QC plots are diagnostic.

---

## Optional utilities (not used by manuscript figures)

The following scripts are bundled for exploratory use. They are not on the
critical path for any manuscript figure or table — the SV-side integration
consumes the merged CSVs from Step 6 directly — but they're useful if you
want to explore the mutation landscape at the gene level.

### `create_mutation_catalog.py`

Aggregates the per-dose merged CSVs from Step 6 into a gene × pattern × dose ×
mutation_type catalog and a wide gene-level summary.

**Input:** directory of `*_dose_*_merged.csv` files
**Output:**
- `mutation_gene_catalog.csv` — long-form: one row per (gene, dose, pattern, mutation_type) combination with counts
- `mutation_gene_summary.csv` — wide: one row per gene with `Mut_0T00`, `Mut_DBS`, `Mut_dA`, etc. breakdown columns

**Command:**
```bash
python create_mutation_catalog.py \
    --mut-dir    merged_analysis \
    --output-dir mutation_catalog
```

Useful for: ranking genes by mutation burden, exploring which genes carry
specific temporal patterns, joining mutation data with external gene lists.

---

## File Format Reference

### 1. seqinfo.txt (from SigProfiler)
```
Sample  Type    Position    Context Ref Alt Strand
d0_W0   SNV     12345      N[C>T]  C   T   0
```

### 2. mutations.pkl (preprocessed)
```
Chromosome, Start, End, Sample, Ref, Alt, Context, Mutation_Strand
Dose, Timepoint, Mutation_Type, Quality_Annotation, MutationID
```

### 3. annotated.pkl (with genomic features)
```
[All columns from mutations.pkl] +
Gene_Name, Feature_Type, Gene_Strand, Distance_To_Gene
```

### 4. mutation_annotations (categorical patterns)
```
MutationID, W0, W1, W2, W3, Category, Pattern
Chromosome, Start, Ref, Alt
```

### 5. merged.csv (complete dataset)
```
[All columns from annotated.pkl] +
[All columns from mutation_annotations] +
Pattern_Category (Radiation-specific/Control-specific/Baseline/Other)
```

---

## Pattern Encoding Reference

### Categorical States (C/T/B/0)
- **C** = Control only
- **T** = Treatment only
- **B** = Both (present in both control and treatment)
- **0** = Absent (neither control nor treatment)

### Example Patterns
- `CBBB` - Present at baseline, persists in both control & treatment
- `0T00` - Treatment-specific at W1 only
- `0TTT` - Treatment-induced, persists across all timepoints
- `C000` - Present at baseline (control), lost at W1

### Pattern Categories
- **Radiation-specific:** `0T00`, `00T0`, `000T`, `0TT0`, `00TT`, `0T0T`, `0TTT`
- **Control-specific:** `C000`, `0C00`, `00C0`, `000C`
- **Baseline:** `CBBB` (present in all timepoints, both conditions)
- **Other patterns:** Any other combination

---

## Quality Control Checkpoints

### After Step 2 (Preprocessing)
```bash
# Check mutation counts
wc -l processed_data/*_mutations.csv

# Check summary statistics
cat summary_data/*_summary.csv
```

### After Step 4 (Annotation)
```bash
# Verify annotation coverage and that real chromosomes survived the join
python -c "
import pandas as pd
df = pd.read_pickle('annotated_mutations/SNV/all_SNV_annotated.pkl')
print(f'Rows: {len(df):,}')
print(f'Annotated (non-Unknown gene): {(df.Gene_Name != \"Unknown\").sum():,}')
print(f'Chromosomes: {df.Chromosome.nunique()} (expect ~24)')
"
```

### After Step 5 (Pattern Assignment)
```bash
# Check pattern distribution AND that chromosomes weren't collapsed to "all"
python -c "
import pandas as pd
df = pd.read_csv('pattern_analysis/SNV/dose_dA/mutation_annotations_dose_dA_SNV.csv')
print(f'Rows: {len(df):,}')
print(f'Chromosome unique: {df.Chromosome.nunique()} (expect ~24, NOT 1)')
print(f'MutationID sample: {df.MutationID.head(3).tolist()}')
print(df.Pattern.value_counts().head(10))
"
```

### After Step 6 (Merging)
```bash
# Verify merge completeness
python -c "
import pandas as pd
df = pd.read_csv('merged_analysis/merged_data/SNV_dose_dA_merged.csv')
print(f'Total: {len(df):,}')
print(f'With genes: {(df.Gene_Name != \"Unknown\").sum():,}')
print(f'With patterns: {(~df.Pattern.isna()).sum():,}')
print(f'Pattern_Group counts:\\n{df.Pattern_Group.value_counts()}')
"
```


## Complete Example Workflow

### Parallelism

**Within a mutation type, Steps 4 → 5 → 6 are a strict chain** — Step 5 needs
Step 4's output and Step 6 needs both. **Across mutation types they're
independent**, so the entire 4→5→6 chain for SNV can run in parallel with the
chain for DBS, ID, and MNS. With 4 cores you get a roughly 4× wall-clock
speedup on the slow steps.

The runtime cost is dominated by Step 4's per-position interval-tree lookup,
and SNV is by far the largest type (~10× DBS) — so the practical limit is RAM,
not CPU. On a machine that can hold all four annotated DataFrames at once,
launch them all; otherwise stagger SNV against the smaller types.


### End-to-end script

```bash
# Step 0: Setup
mkdir -p {vcf_files,sigprofiler_output,processed_data,annotations,annotated_mutations,pattern_analysis,merged_analysis,sankey_output}

# Step 1: SigProfiler (if starting from VCF)
python sigprofiler.py -i vcf_files -o sigprofiler_output -r GRCh38

# Step 2: Preprocess mutations (writes processed_data/all_<TYPE>_mutations.pkl)
python mutation_preprocessing.py \
    --input-dir vcf_files/output/vcf_files \
    --chromosome all \
    --output processed_data \
    --summary summary_data

# Step 3: Prepare annotations (one-time, shared across all types)
python annotation_preprocessing.py --build hg38 --annotation-dir annotations

# Steps 4-6: per-type chain. Loop sequentially or background-run for parallelism.
for mut_type in SNV DBS ID MNS; do
    # Step 4: annotate
    python mutation_annotation.py \
        -m processed_data/all_${mut_type}_mutations.pkl \
        -a annotations \
        -b hg38 \
        -o annotated_mutations/${mut_type}

    # Rename Step 4's fixed output filename to what Step 5 expects to glob
    mv annotated_mutations/${mut_type}/annotated_mutations.pkl \
       annotated_mutations/${mut_type}/all_${mut_type}_annotated.pkl
    mv annotated_mutations/${mut_type}/annotated_mutations.csv \
       annotated_mutations/${mut_type}/all_${mut_type}_annotated.csv

    # Step 5: assign temporal patterns (one categorical + one binary file per dose)
    python mutation_pattern_assignment.py \
        -i annotated_mutations/${mut_type} \
        -o pattern_analysis/${mut_type} \
        -m ${mut_type}

    # Step 6: merge annotations with patterns. Dose names match dirs from Step 5.
    python merge_annotation.py \
        --annotated-dir annotated_mutations \
        --pattern-dir pattern_analysis \
        --output-dir merged_analysis \
        --mutation-types ${mut_type} \
        --doses dA dB dC dD dE
done

# Step 7: Sankey diagrams. Compute trajectories once per type (writes both
# per-dose and combined JSONs), then render the combined manuscript figure
# from each type's combined trajectories JSON.
for mut_type in SNV DBS ID MNS; do
    python compute_sankey.py \
        --input      annotated_mutations/${mut_type}/all_${mut_type}_annotated.pkl \
        --output-dir sankey_flows/${mut_type}

    python render_sankey.py \
        --trajectories-json sankey_flows/${mut_type}/combined/all_chromosomes_trajectories.json \
        --output            sankey_figures/${mut_type}_combined.png \
        --title             "Temporal Dynamics for ${mut_type}" \
        --subtitle          "All Doses, All Chromosomes"
done

# Step 8: Downstream analyses (use the merged CSVs)
python your_analysis_script.py -i merged_analysis/merged_data/SNV_dose_dD_merged.csv
```

---

## Required Python Packages

All Python dependencies are pinned in `requirements.txt`:

```bash
pip install -r requirements.txt
```

For GPU acceleration of `sigprofiler.py`, also install a CUDA-matched PyTorch
build (see **Environment Setup** at the top of this README).

---

## Summary

This pipeline takes you from raw VCF files to comprehensive mutation analysis:

1. ✅ **VCF → seqinfo.txt** (SigProfiler)
2. ✅ **seqinfo.txt → structured mutations** (preprocessing)
3. ✅ **Genomic annotations** (one-time setup)
4. ✅ **Mutation + genomic features** (annotation)
5. ✅ **Temporal patterns** (pattern assignment)
6. ✅ **Complete merged dataset** (merge)
7. ✅ **Visualizations** (Sankey diagrams)
8. ✅ **Ready for analysis** (your custom scripts)

Each step is modular, reusable, and produces both machine-readable (`.pkl`) and human-readable (`.csv`) outputs.
