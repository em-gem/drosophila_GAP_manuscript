#!/bin/bash
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    echo "Example: $0 /path/to/fastq/directory"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="results/dropseq_processed"
mkdir -p ${OUTPUT_DIR}

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist"
    exit 1
fi

# Check if input files exist with expected pattern
if ! ls ${INPUT_DIR}/*_R1.fastq.gz >/dev/null 2>&1; then
    echo "Error: No files matching *_R1.fastq.gz found in $INPUT_DIR"
    exit 1
fi

# Process each sample (assuming files are named sample_R1.fastq.gz and sample_R2.fastq.gz)
for r1_file in ${INPUT_DIR}/*_R1.fastq.gz; do
    sample=$(basename $r1_file _R1.fastq.gz)
    r2_file="${INPUT_DIR}/${sample}_R2.fastq.gz"
    
    echo "Processing sample ${sample}..."
        
    # 1. Convert FASTQ to uBAM
    picard FastqToSam \
        F1=${INPUT_DIR}/${sample}_R1.fastq.gz \
        F2=${INPUT_DIR}/${sample}_R2.fastq.gz \
        O=${OUTPUT_DIR}/${sample}_unaligned.bam \
        SAMPLE_NAME=${sample}

    # 2. Tag cell barcodes
    dropseq TagBamWithReadSequenceExtended \
        INPUT=${OUTPUT_DIR}/${sample}_unaligned.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_tagged_Cell.bam \
        SUMMARY=${OUTPUT_DIR}/${sample}_tagged_Cellular.bam_summary.txt \
        BASE_RANGE=1-12 \
        BASE_QUALITY=10 \
        BARCODED_READ=1 \
        DISCARD_READ=False \
        TAG_NAME=XC \
        NUM_BASES_BELOW_QUALITY=1

    rm ${OUTPUT_DIR}/${sample}_unaligned.bam  

    # 3. Tag molecular barcodes
    dropseq TagBamWithReadSequenceExtended \
        INPUT=${OUTPUT_DIR}/${sample}_tagged_Cell.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_tagged_CellMolecular.bam \
        SUMMARY=${OUTPUT_DIR}/${sample}_tagged_Molecular.bam_summary.txt \
        BASE_RANGE=13-20 \
        BASE_QUALITY=10 \
        BARCODED_READ=1 \
        DISCARD_READ=True \
        TAG_NAME=XM \
        NUM_BASES_BELOW_QUALITY=1

    rm ${OUTPUT_DIR}/${sample}_tagged_Cell.bam 

    # 4. Filter BAM
    dropseq FilterBam \
        TAG_REJECT=XQ \
        INPUT=${OUTPUT_DIR}/${sample}_tagged_CellMolecular.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_tagged_filtered.bam

    rm ${OUTPUT_DIR}/${sample}_tagged_CellMolecular.bam 

    # 5. Trim SMART adapter
    dropseq TrimStartingSequence \
        INPUT=${OUTPUT_DIR}/${sample}_tagged_filtered.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_trimmed_smart.bam \
        OUTPUT_SUMMARY=${OUTPUT_DIR}/${sample}_adapter_trimming_report.txt \
        SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
        MISMATCHES=0 \
        NUM_BASES=5

    rm ${OUTPUT_DIR}/${sample}_tagged_filtered.bam

    # 6. Trim polyA
    dropseq PolyATrimmer \
        INPUT=${OUTPUT_DIR}/${sample}_trimmed_smart.bam \
        OUTPUT=${OUTPUT_DIR}/${sample}_polyA_filtered.bam \
        OUTPUT_SUMMARY=${OUTPUT_DIR}/${sample}_polyA_trimming_report.txt \
        MISMATCHES=0 \
        NUM_BASES=6 \
        USE_NEW_TRIMMER=true

   rm ${OUTPUT_DIR}/${sample}_trimmed_smart.bam 

    # 7. Convert to FASTQ for STAR
    picard SamToFastq \
        INPUT=${OUTPUT_DIR}/${sample}_polyA_filtered.bam \
        FASTQ=${OUTPUT_DIR}/${sample}_for_star.fastq

    echo "Finished processing ${sample}"
done
