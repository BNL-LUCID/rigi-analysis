#!/bin/bash
# Script to create the Conda environment for AnnotSV
# Kept in its own environment due to TCL/openblas dependencies

set -e

echo "Creating 'annotsv' Conda environment..."
conda create -y -n annotsv -c bioconda annotsv
conda run -n annotsv conda install -y -c conda-forge openblas

echo ""
echo "AnnotSV environment created successfully!"
echo "To activate the environment, run:"
echo "conda activate annotsv"
echo ""
echo "Note: You must also download the AnnotSV database (GRCh38) from https://lbgi.fr/AnnotSV/ and unpack it."
echo "Provide the path to the unpacked database to the pipeline scripts using the --annotationsDir flag or ANNOTSV_DB variable."
