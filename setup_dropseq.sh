#!/bin/bash
set -e

# Clone Drop-seq tools
git clone https://github.com/broadinstitute/Drop-seq.git
cd Drop-seq

# Build the executable jarfile
./gradlew installDist

# Create symbolic link to make dropseq available in PATH
ln -s "$(pwd)/dropseq/build/install/dropseq/bin/dropseq" "$CONDA_PREFIX/bin/dropseq"

# Test installation
dropseq --help
