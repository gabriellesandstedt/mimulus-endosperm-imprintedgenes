#!/bin/bash

# Path to the GFF3 file
gff_file="path/to/your/gff3/file.gff3"

# Path to the output file
output_file="gene_lengths.txt"

# Extract gene features from the GFF3 file
awk '$3 == "gene" { print $1"\t"$4"\t"$5"\t"$9 }' MguttatusvarIM62v3.1.primaryTrs.gff3 > genes.bed

awk -v OFS='\t' '{split($4,a,"="); split(a[2],b,";"); print $1, $2, $3, b[1]}' genes.bed > genes2.bed

awk -v OFS='\t' '{print $1, $3 - $2, $4}' genes2.bed > gene_lengths.txt



awk -v OFS='\t' '{print $1, $2}' gene_lengths_unique_shared_til_MEGs.txt > gene_lengths_unique_shared_til_MEGs2.txt
