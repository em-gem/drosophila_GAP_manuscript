#!/bin/bash
set -e

# Check if references directory exists
if [ ! -d "references" ]; then
    echo "Error: references directory not found. Run download_reference.sh first."
    exit 1
fi

echo "Creating sequence dictionary..."
picard CreateSequenceDictionary \
  R=references/dmel-all-chromosome-r6.63.with_transgenes.fasta \
  O=references/dmel-all-chromosome-r6.63.with_transgenes.dict

echo "Converting GTF to refFlat..."
dropseq ConvertToRefFlat \
  ANNOTATIONS_FILE=references/dmel-all-r6.63.with_transgenes.gtf \
  OUTPUT=references/dmel-all-r6.63.with_transgenes.refFlat \
  SEQUENCE_DICTIONARY=references/dmel-all-chromosome-r6.63.with_transgenes.dict

echo "Reducing GTF..."
dropseq ReduceGtf \
  GTF=references/dmel-all-r6.63.with_transgenes.gtf \
  OUTPUT=references/dmel-all-r6.63.reduced.with_transgenes.gtf \
  SEQUENCE_DICTIONARY=references/dmel-all-chromosome-r6.63.with_transgenes.dict \
  ENHANCE_GTF=true

echo "Creating interval files..."
mkdir -p references/intervals

dropseq CreateIntervalsFiles \
  SEQUENCE_DICTIONARY=references/dmel-all-chromosome-r6.63.with_transgenes.dict \
  REDUCED_GTF=references/dmel-all-r6.63.reduced.with_transgenes.gtf \
  OUTPUT=references/intervals \
  PREFIX=dmel-all-r6.63 \
  MT_SEQUENCE=mitochondrion_genome

echo "Reference preparation complete!"
