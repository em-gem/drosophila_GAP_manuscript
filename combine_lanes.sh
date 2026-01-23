#!/bin/bash
set -e

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    echo "Example: $0 /Volumes/Data/RNASeq/allFastq data/merged_fastq"
    exit 1
fi

# Define directories from arguments
FASTQ_DIR="$1"
OUTPUT_DIR="$2"

# Check if input directory exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: Input directory $FASTQ_DIR does not exist"
    exit 1
fi

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Find unique sample names (F1, F2, d1)
for sample in F1 F2 d1; do
    echo "Processing sample ${sample}..."
    
    # Merge R1 files (cell barcodes)
    cat ${FASTQ_DIR}/dropseq-${sample}_S*_L00{1,2,3,4}_R1_001.fastq.gz > ${OUTPUT_DIR}/${sample}_R1.fastq.gz
    
    # Merge R2 files (transcripts)
    cat ${FASTQ_DIR}/dropseq-${sample}_S*_L00{1,2,3,4}_R2_001.fastq.gz > ${OUTPUT_DIR}/${sample}_R2.fastq.gz
    
    echo "Finished merging lanes for sample ${sample}"
done
