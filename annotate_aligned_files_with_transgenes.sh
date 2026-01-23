#!/bin/bash
set -e

if [ $# -ne 3 ]; then
    echo "Usage: $0 <aligned_bam_dir> <tagged_bam_dir> <output_dir>"
    echo "Example: $0 results/star_aligned_with_transgenes results/dropseq_processed results/final_dge"
    exit 1
fi

ALIGNED_DIR="$1"
TAGGED_DIR="$2"
OUTPUT_DIR="$3"
TMP_DIR="./tmp"

mkdir -p ${OUTPUT_DIR}
mkdir -p ${TMP_DIR}  # make sure TMP_DIR exists


for aligned_bam in ${ALIGNED_DIR}/*_Aligned.sortedByCoord.out.bam; do
    sample=$(basename $aligned_bam _Aligned.sortedByCoord.out.bam)
    echo "Processing ${sample}..."

    # 1. Sort aligned BAM by queryname
    if [ ! -f ${OUTPUT_DIR}/${sample}_queryname_sorted.bam ]; then
        picard SortSam \
            I=${aligned_bam} \
            O=${OUTPUT_DIR}/${sample}_queryname_sorted.bam \
            SORT_ORDER=queryname \
            TMP_DIR=${TMP_DIR}
    else
        echo "Skipping SortSam for ${sample}, output already exists."
    fi


    # 2. Merge aligned BAM with cell/molecular barcodes
    if [ ! -f ${OUTPUT_DIR}/${sample}_merged.bam ]; then
        picard MergeBamAlignment \
            REFERENCE_SEQUENCE=references/dmel-all-chromosome-r6.63.with_transgenes.fasta \
            UNMAPPED_BAM=${TAGGED_DIR}/${sample}_polyA_filtered.bam \
            ALIGNED_BAM=${OUTPUT_DIR}/${sample}_queryname_sorted.bam \
            OUTPUT=${OUTPUT_DIR}/${sample}_merged.bam \
            INCLUDE_SECONDARY_ALIGNMENTS=false \
            PAIRED_RUN=false \
            TMP_DIR=${TMP_DIR}
    else
        echo "Skipping MergeBamAlignment for ${sample}, output already exists."
    fi

    # 3. Tag with gene exons
    if [ ! -f ${OUTPUT_DIR}/${sample}_gene_tagged.bam ]; then
        dropseq TagReadWithGeneFunction \
             I=${OUTPUT_DIR}/${sample}_merged.bam \
             O=${OUTPUT_DIR}/${sample}_gene_tagged.bam \
             ANNOTATIONS_FILE=references/dmel-all-r6.63.with_transgenes.gtf \
             TMP_DIR=${TMP_DIR}
    else
	echo "Skipping TagReadWithGeneExonFunction for ${sample}, output already exists."
    fi

    # 4. Create Digital Gene Expression matrix
    dropseq DigitalExpression \
         I=${OUTPUT_DIR}/${sample}_gene_tagged.bam \
         O=${OUTPUT_DIR}/${sample}_100k_barcodes_dge.txt.gz \
         CELL_BARCODE_TAG=XC \
         SUMMARY=${OUTPUT_DIR}/${sample}_100k_barcodes_dge_summary.txt \
         NUM_CORE_BARCODES=100000 \
         TMP_DIR=${TMP_DIR}

    # Clean up intermediate files
    #rm ${OUTPUT_DIR}/${sample}_queryname_sorted.bam
    #rm ${OUTPUT_DIR}/${sample}_merged.bam

    echo "Finished processing ${sample}"
done

