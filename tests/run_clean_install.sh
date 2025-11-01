#!/bin/bash

# Script to test clean installation of genePackRat
# This ensures proper environment isolation on HPC systems

echo "Setting up clean test environment for genePackRat..."

# Clear any existing R environment variables
unset R_LIBS_USER
unset R_LIBS_SITE
unset R_LIBS

# Ensure we're in the correct directory
cd "$(dirname "$0")"

# Check if conda environment exists
if ! conda env list | grep -q "test_genePackRat"; then
    echo "Creating conda environment from clean_install_env.yml..."
    conda env create -f clean_install_env.yml
fi

# Activate the environment and run the test
echo "Activating test_genePackRat environment..."
eval "$(conda shell.bash hook)"
conda activate test_genePackRat

# Export conda library path for R
export R_LIBS="${CONDA_PREFIX}/lib/R/library"

echo "Running installation test..."
Rscript clean_install.R

echo "Test complete!"