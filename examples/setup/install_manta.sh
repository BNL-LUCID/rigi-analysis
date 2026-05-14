#!/bin/bash
# Script to create the Conda environment for Manta SV calling
# Manta 1.6.0 requires Python 2.7, so it is kept in a separate environment.

set -e

echo "Creating 'sv_tools' Conda environment for Manta..."
conda create -y -n sv_tools -c conda-forge -c bioconda python=2.7 manta samtools

echo "Verifying Manta installation..."
conda run -n sv_tools manta --version

echo ""
echo "Installation complete! To activate the environment, run:"
echo "conda activate sv_tools"
