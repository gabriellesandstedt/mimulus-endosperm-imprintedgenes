#!/bin/bash

# path to bed file, that contains 10 columns. col 1 | chr, col 2| position start, col 3| position end, col 4 | gene name... col 8 | genome feature
# bed file was generated from gff/gtf 
bed_file="genes.bed"

# extract gene features and calculate lengths
awk '$8 ~ /gene/' "$bed_file" | awk -F'\t' '{print $1, $3 - $2, $4}' > gene_lengths.txt
