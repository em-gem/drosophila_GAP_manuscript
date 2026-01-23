#!/bin/bash
set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    echo "Example: $0 results/dropseq_processed results/star_aligned_with_transgenes"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
mkdir -p ${OUTPUT_DIR}

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist"
    exit 1
fi

# File path handling is pretty bad here..
if [ ! -d "references/star_index_with_transgenes" ]; then
    echo "Error: STAR index not found in references/star_index_with_transgenes"
    exit 1
fi

for fastq in ${INPUT_DIR}/*_for_star.fastq; do
    sample=$(basename $fastq _for_star.fastq)
    echo "Aligning sample ${sample}..."
    
    STAR --genomeDir references/star_index_with_transgenes \
        --readFilesIn ${fastq} \
        --outFileNamePrefix ${OUTPUT_DIR}/${sample}_ \
        --runThreadN 8 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1

    echo "Finished aligning ${sample}"
done
