# Mutation Analysis Pipeline

Point-mutation half of the LUCID radiation-mutagenesis pipeline. Takes 
VCF files through SigProfiler matrix generation, genomic annotation, 
and temporal pattern assignment, producing the per-dose merged tables 
that the SV-side scripts consume.

## Pipeline
VCF files
↓ [1] sigprofiler.py            → seqinfo.txt
↓ [2] mutation_preprocessing.py → {chr}{type}mutations.pkl
↓ [3] annotation_preprocessing.py → hg38 interval trees (one-time)
↓ [4] mutation_annotation.py    → all{type}annotated.pkl
↓ [5] mutation_pattern_assignment.py → mutation_annotations_dose.csv
↓ [6] merge_annotation.py       → {type}dose{dose}_merged.csv  ← consumed by SV side
↓ [7] compute_sankey.py + render_sankey.py → Fig 2, S1

Steps 4 → 5 → 6 are a strict chain within a mutation type; across types 
they're independent, so SNV/DBS/ID/MNS can run concurrently if RAM allows.

## Setup

```bash
conda create -n lucid_mut python=3.9 -y
conda activate lucid_mut
pip install -r requirements.txt
```

For offline genome installation (HPC/SciServer environments), pre-download 
GRCh38.tar.gz from the AlexandrovLab FTP and pass `--offline-genome` to 
step 1.

GPU acceleration of step 1 is optional; CUDA-matched PyTorch must be 
installed separately if used.

## End-to-end example

```bash
# One-time setup
python annotation_preprocessing.py --build hg38 --annotation-dir annotations

# Per-VCF run
python sigprofiler.py -i vcf_files -o sigprofiler_output -r GRCh38

python mutation_preprocessing.py \
    --input-dir vcf_files/output/vcf_files \
    --output processed_data \
    --summary summary_data

# Steps 4-6 per mutation type
for mut_type in SNV DBS ID MNS; do
    python mutation_annotation.py \
        -m processed_data/all_${mut_type}_mutations.pkl \
        -a annotations -b hg38 \
        -o annotated_mutations/${mut_type}
    
    # Step 5 expects all_<TYPE>_annotated.pkl naming
    mv annotated_mutations/${mut_type}/annotated_mutations.pkl \
       annotated_mutations/${mut_type}/all_${mut_type}_annotated.pkl
    
    python mutation_pattern_assignment.py \
        -i annotated_mutations/${mut_type} \
        -o pattern_analysis/${mut_type} \
        -m ${mut_type}
    
    python merge_annotation.py \
        --annotated-dir annotated_mutations \
        --pattern-dir pattern_analysis \
        --output-dir merged_analysis \
        --mutation-types ${mut_type} \
        --doses dA dB dC dD dE
done

# Step 7: Sankey figures
for mut_type in SNV DBS ID MNS; do
    python compute_sankey.py \
        --input annotated_mutations/${mut_type}/all_${mut_type}_annotated.pkl \
        --output-dir sankey_flows/${mut_type}
    
    python render_sankey.py \
        --trajectories-json sankey_flows/${mut_type}/combined/all_chromosomes_trajectories.json \
        --output sankey_figures/${mut_type}_combined.png \
        --title "Temporal Dynamics for ${mut_type}" \
        --subtitle "All Doses, All Chromosomes"
done
```

## Per-script documentation

Each script's `--help` and module docstring describe its arguments, 
inputs, outputs, and edge cases. Highlights:

- **`mutation_preprocessing.py`**: searches `<input-dir>/{SNV,DBS,ID,MNS}/` 
  recursively for `*_seqinfo.txt`. SigProfiler's nested layout works directly.
- **`mutation_annotation.py`**: produces single fixed-name output 
  (`annotated_mutations.pkl`); rename to `all_<TYPE>_annotated.pkl` before 
  step 5 globs for it.
- **`compute_sankey.py`**: emits both trajectory JSONs (preferred) and 
  legacy pairwise flow JSONs. The combined view is computed directly from 
  the full dataset — do not sum per-dose JSONs (controls would be counted 5×).

## Parallel helpers

`parallel_annotate.py` and `parallel_pattern_assignment.py` wrap steps 4 
and 5 with chunk-level and dose-level parallelism respectively. See their 
`--help` and `docs/PARALLEL.md` for caveats around RAM scaling and per-chunk 
QC outputs.

## Reference docs

- [`docs/PATTERNS.md`](docs/PATTERNS.md) — pattern encoding and category definitions
- [`docs/FILE_FORMATS.md`](docs/FILE_FORMATS.md) — column schemas at each step
- [`docs/QC.md`](docs/QC.md) — sanity-check snippets after each step
- [`docs/PARALLEL.md`](docs/PARALLEL.md) — parallel helper caveats
