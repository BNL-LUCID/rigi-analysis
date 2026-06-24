# LUCID — Spatio-Temporal Analysis of Radiation-Induced Genomic Instability

Analytical code accompanying:

> **A Spatio-Temporal Analysis Framework for Characterizing
> Radiation-Induced Genomic Instability**
> Chopra K., Cucinell C., Weinberg R., Forrester S., Brettin T.,
> Kilic O. O., Yoon B. J. (2026).

The pipeline analyses whole-genome sequencing of HUVEC cells exposed
to chronic low-dose gamma radiation (0.20–2.62 mGy/hr) over three
weeks. It produces the temporal pattern catalogues, breakpoint-proximal
mutation enrichment statistics, dose-stratified analyses, and figures
presented in the paper.

## Repository Layout & Pipelines

The analysis framework is divided into two primary pipelines that meet at an intermediate artefact (dose-merged mutation tables):

1. **Mutation Analysis Pipeline**
   - **Location:** `src/rigi_analysis/mutation/`
   - **Flow:** `Mutect2` → `SigProfiler` → temporal patterns → `Sankey` visualizations.
   - **Outputs:** `*_dose_*_merged.csv` tables which are consumed by the SV pipeline.
   - **Documentation:** [Mutation Pipeline Instructions](mutation/README.md)

2. **Structural-Variant (SV) Pipeline**
   - **Location:** `src/rigi_analysis/sv/`
   - **Flow:** `Manta` → `AnnotSV` → temporal patterns → breakpoint enrichment → dose stratification.
   - **Inputs:** Reads the merged CSVs produced from the mutation half.
   - **Documentation:** [SV Pipeline Instructions](sv/README.md)

**Execution Order:** You must run the `mutation` pipeline first to generate the
required intermediate tables, followed by the `sv` pipeline. Each documentation
file provides specific setup instructions, environmental requirements, and
end-to-end examples.

---

## Workflow Entry Points

Three workflow entry points are provided as console scripts. Each launches async pipeline 
stages via [RADICAL](https://radical-cybertools.github.io) tools, backed by the execution 
engine configured in [`backend.py`](../src/rigi_analysis/workflows/backend.py) 
(default: `DaskExecutionBackend`).

| Command | Module | Description |
|---------|--------|-------------|
| `rigi-analysis-workflow` | `full_workflow:main` | Runs `mutation` → `sv` pipelines sequentially in a single session |
| `rigi-analysis-workflow-mutation` | `mutation_workflow:main` | Runs the `mutation` pipeline only |
| `rigi-analysis-workflow-sv` | `sv_workflow:main` | Runs the `sv` pipeline only |

All three accept the same CLI flags:

```
-c, --config-file   JSON configuration file (required)
-o, --output-dir    Top-level output directory
```

The per-stage script dispatcher is available via `rigi-analysis-run <script_name>` for manual step-by-step execution.

---

## Workflow Data Structures

### Configuration JSON

All three workflows read a single JSON configuration file (`-c`). All pipeline parameters are **top-level keys**. The only special nested key is `run_description`, which is passed to the execution backend:

```json
{
  "annotations_dir": "./annotations",             // required
  "vcf_dir": "./vcf_files",                       // required
  "genome_build": "hg38",
  "sigprofiler_ref": "GRCh38",
  "mutation_types": ["SNV", "DBS", "ID", "MNS"],  // required
  "doses": ["dA", "dB", "dC", "dD", "dE"],        // required
  
  "annotsv_dir": "./annotsv",                     // required
  "annotsv_control_dir": "./annotsv_control",     // required
  "sv_tolerance": 1000,
  "windows": "10,25,50,100",
  "single_window": 10,
  "mega_threshold": 50000000,
  "gnomad_file": "./gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz",
  "run_description": {}
}
```

#### Mutation Workflow Keys

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `annotations_dir` | `str` | `"<cwd>/annotations"` | Directory for interval-tree annotation caches |
| `vcf_dir` | `str` | `"<cwd>/vcf_files"` | Input directory containing raw VCF files (Mutect2 output) |
| `genome_build` | `str` | `"hg38"` | Genome build for annotation preprocessing |
| `sigprofiler_ref` | `str` | `"GRCh38"` | SigProfiler reference genome identifier |
| `mutation_types` | `list[str]` | `["SNV","DBS","ID","MNS"]` | Mutation types to process |
| `doses` | `list[str]` | `["dA","dB","dC","dD","dE"]` | Dose labels for merge step |

#### SV Workflow Keys

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `annotsv_dir` | `str` | `"<cwd>/annotsv"` | Input directory of AnnotSV-annotated TSVs (radiation) |
| `annotsv_control_dir` | `str` | `"<cwd>/annotsv_control"` | Input directory of AnnotSV-annotated TSVs (control) |
| `sv_tolerance` | `int` | `1000` | Breakpoint position tolerance (bp) for cross-timepoint SV matching |
| `windows` | `str` | `"10,25,50,100"` | Comma-separated breakpoint-proximal window sizes (bp) |
| `mutation_merged_dir` | `str` | `"../out_mutation/merged_data"` | Merged mutation CSVs from mutation pipeline (not needed for a full workflow run) |
| `gnomad_file` | `str` | `"<cwd>/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"` | gnomAD pLI constraint file |
| `single_window` | `int` | `10` | Single window size for dose-stratified and concordance analyses |
| `mega_threshold` | `int` | `50000000` | Inversion size threshold (bp) for mega-inversion classification |

#### Full Workflow

Uses both mutation and SV keys from the same config. The full workflow automatically wires `mutation_merged_dir` in the SV config to the mutation pipeline's merged output — no manual path threading needed.

#### `run_description` — Backend Configuration

Nested dict passed directly to the execution backend. Contents are backend-specific (e.g., resource, runtime, n_workers, threads_per_worker, etc.).

---

## Mutation Workflow — I/O Directory Structure

**Command:** `rigi-analysis-workflow-mutation -c config.json -o out_mutation`

```
annotations/                         # Interval-tree caches (stage 1a)
vcf_files/                           # Raw VCFs from Mutect2 pipeline 
└── output/vcf_files/                # Parsed SigProfiler outputs into per-type pickles

out_mutation/
├── sigprofiler_output/              # SigProfiler matrices & seqinfo (stage 1b)
├── summary_data/                    # Per-type annotation summaries (stage 2)
├── processed_data/                  # Parsed per-type pickles (stage 2)
│   ├── all_<TYPE>_mutations.pkl
│   └── ...
├── annotated_mutations/             # Genomic-annotated pickles (stage 3)
│   ├── <TYPE>/all_<TYPE>_annotated.pkl
│   └── ...
├── pattern_analysis/                # Temporal pattern CSVs (stage 4)
│   ├── <TYPE>/mutation_annotations_dose_<DOSE>_<TYPE>.csv
│   └── ...
├── merged_analysis/                 # ★ Handover to SV pipeline (stage 5)
│   ├── <TYPE>_dose_<DOSE>_merged.csv
│   └── ...
├── sankey_flows/                    # Trajectory JSONs (stage 6)
│   ├── <TYPE>/combined/all_chromosomes_trajectories.json
│   └── ...
└── sankey_figures/                  # Sankey PNGs (stage 7)
    ├── <TYPE>_combined.png
    └── ...
```

### DAG

```
annotation_preprocessing ──┐
sigprofiler ──► preprocessing ──┐
                                └──► annotation(×N) ──┬──► pattern ──► merge
                                                      └──► compute_sankey ──► render_sankey
```

Stages 1a (`annotation_preprocessing`) and 1b (`sigprofiler`) run in parallel. Stages 3–7 are instantiated per mutation type and run concurrently. Within each type, the Sankey branch (6→7) runs in parallel with the pattern→merge branch (4→5).

---

## SV Workflow — I/O Directory Structure

**Command:** `rigi-analysis-workflow-sv -c config.json -o out_sv`

**Required inputs (from mutation pipeline):** `mutation_merged_dir` pointing to the merged CSVs above (`<mutation_dir>/merged_data`, and specifically `DBS` subdirectory for INV-DBS enrichment analyses).

```
annotsv_passed/                       # PASS-filtered AnnotSV (step 1a)
annotsv_passed_control/               # PASS-filtered control AnnotSV (step 1b)

out_sv/
├── sv_temporal/                      # Temporal catalog (step 2)
│   └── sv_temporal_catalog.csv       #   → consumed by steps 5, 6
├── figure3_sv_landscape.png          # SV landscape figure (step 3)
├── figure_s2_repeat_analysis.png     # Repeat analysis figure (step 4)
├── sv_correlation/                   # Breakpoint enrichment results (step 5)
│   └── DBS_overall/enrichment_by_window.csv
├── sv_type_decay/                    # SV-type decay results (step 6)
│   └── DBS/sv_type_enrichment_decay.csv
├── inv_size/                         # INV-size analysis (step 7)
│   ├── inv_size.png
│   ├── size_class_results.csv
│   └── mega_test_results.csv
├── figure4_inv_dbs_unified.png       # INV-DBS coupling figure (step 8)
├── figure_temporal_dynamics.png      # Temporal concordance figure (step 10)
├── dose_stratified/                  # Dose-stratified pairs (step 11)
│   ├── genes_high_dose.csv
│   ├── genes_low_dose.csv
│   ├── inv_dbs_pairs_low.csv
│   └── inv_dbs_pairs_high.csv
├── annotated_genes.csv               # MyGene.info + gnomAD annotations (step 12)
├── categorized_genes.csv             # 6 functional categories (step 13)
└── figure5_dose_response/            # Dose-response figure (step 14)
```

### DAG

```
filter_pass_rad ──┬──► sv_temporal ──┬──► correlation ──┐
                  │                  └──► type_decay ───┤
                  ├──► sv_landscape (terminal)          ├──► sv_mut_viz (Fig 4)
                  ├──► inv_size ────────────────────────┘
                  ├──► concordance ──► concordance_viz (Fig S4)
                  └──► dose_stratified ──► fetch ──► categorise ──► dose_viz (Fig 5)
filter_pass_ctl ──┴──► repeat_analysis (terminal)
```

Steps 1a and 1b run in parallel. Steps 3, 4, 7, 9, 11 all depend only on step 1 and launch concurrently. Steps 5–6 require the temporal catalog from step 2. Step 8 waits for steps 5, 6, and 7.

---

## Full Workflow — I/O Directory Structure

**Command:** `rigi-analysis-workflow -c config.json -o out_full`

Composes both pipelines sequentially. The `mutation` pipeline runs to completion first, then its `merged_data/` directory is automatically wired through `mutation_merged_dir` for the `sv` pipeline.

```
out_full/
├── out_mutation/        # Full Mutation directory tree (see above)
│   └── merged_data/     # ★ auto-wired to SV pipeline
└── out_sv/              # Full SV directory tree (see above)
```

---

## Setup

A unified Conda environment handles all Python dependencies for the `rigi_analysis` package.

```bash
# From the repository root
conda env create -f environment.yml
conda activate rigi_analysis
pip install -e ".[dev]"
```

For GPU-accelerated SigProfiler signature extraction, install the optional extra:

```bash
pip install -e ".[gpu]"   # adds torch, torchvision, torchaudio
```

### Python Dependencies

Key dependencies (all installed automatically):

| Group | Packages |
|-------|----------|
| Core scientific stack | `numpy>=1.23`, `pandas>=1.5`, `scipy>=1.10` |
| Parallelism | `dask[dataframe]>=2023.1.0`, `joblib` |
| Genomic lookups / REST | `intervaltree>=3.1.0`, `requests>=2.28` |
| Plotting | `matplotlib>=3.6`, `seaborn>=0.12` |
| Mutational signatures | `SigProfilerMatrixGenerator>=1.2.25`, `SigProfilerExtractor>=1.1.24` |
| Workflow orchestration | `radical.asyncflow>=0.3.1`, `rhapsody-py[dask]>=0.2.0` |
| GPU (optional) | `torch`, `torchvision`, `torchaudio` (via `.[gpu]` extra) |

### Third-Party Tools & Data

Several external bioinformatic tools and reference files are required upstream:
- **Manta 1.6.0**: Somatic SV calling. Requires `samtools>=1.10` and its own Python 2.7 environment.
- **AnnotSV 3.3**: SV annotation (must be run in `full` mode).
- **Mutect2**: Somatic point mutation calling.
- **Reference Data**:
  - gnomAD v2.1.1 LoF metrics (`gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz`).
  - SigProfiler reference genomes (e.g., GRCh38).

*Setup scripts for Manta and AnnotSV Conda environments are in `examples/setup/`. See [sv/install.md](sv/install.md) for details.*

## Execution Backend

The workflow execution backend is configured centrally in [`backend.py`](../src/rigi_analysis/workflows/backend.py). Change the `_BACKEND_CLASS` assignment to switch all workflows at once:

| Backend | Import | Use case |
|---------|--------|----------|
| `DaskExecutionBackend` | `rhapsody.backends` | **Default.** Dask Distributed cluster |
| `RadicalExecutionBackend` | `rhapsody.backends` | RADICAL-Pilot on HPC (requires `rhapsody-py[radical_pilot]`) |
| `DragonExecutionBackendV3` | `rhapsody.backends` | DragonHPC runtime |
| `Concurrent` | `rhapsody.backends` | Thread / process pool |

## Citation

If you use this pipeline, please cite the manuscript above. Upstream tools (Mutect2, SigProfiler, Manta, AnnotSV) should also be cited separately — see the manuscript References section.

## Funding

This work was supported by the U.S. Department of Energy, Office of Science, Biological and Environmental Research program (B&R# KP1601017, FWP# CC140).

## Out of Scope

Raw sequencing data, alignment, and upstream variant calling are not included. See the manuscript Methods (§4) and the BioProject accession in Data Access for those.
