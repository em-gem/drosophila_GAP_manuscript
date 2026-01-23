#!/bin/bash
set -e

mkdir -p references/star_index_with_transgenes

STAR --runMode genomeGenerate \
    --genomeDir references/star_index_with_transgenes \
    --genomeFastaFiles references/dmel-all-chromosome-r6.63.with_transgenes.fasta \
    --sjdbGTFfile references/dmel-all-r6.63.with_transgenes.gtf \
    --runThreadN 8 \
    --genomeSAindexNbases 12

