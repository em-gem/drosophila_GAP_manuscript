#!/bin/bash
set -e

mkdir -p references

wget --no-check-certificate -P references https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.63_FB2025_02/fasta/dmel-all-chromosome-r6.63.fasta.gz
wget --no-check-certificate -P references https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.63_FB2025_02/gtf/dmel-all-r6.63.gtf.gz
gunzip references/*.gz

