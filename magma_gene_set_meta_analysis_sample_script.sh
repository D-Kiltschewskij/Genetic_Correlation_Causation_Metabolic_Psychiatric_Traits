#!/bin/bash

cd </path/to/working/directory>

# Identify genes with data for both the metabolite and psychiatric condition (using output from gene-level meta-analysis)
awk -F " " 'NR>1{print $1}' GCST90451108_ADHD_PGC_2022_meta_common_genes.txt > GCST90451108_ADHD_PGC_2022_meta_common_genes_ids.txt

# Make new gene annotation files for each trait, subset to common genes
head -2 GCST90451108.genes.annot > GCST90451108_GCST90451108_ADHD_PGC_2022_common_genes.genes.annot
head -2 ADHD_PGC_2022.genes.annot > ADHD_PGC_2022_GCST90451108_ADHD_PGC_2022_common_genes.genes.annot

while read -r gene
do
 grep -P "^${gene}\s" GCST90451108.genes.annot >> GCST90451108_GCST90451108_ADHD_PGC_2022_common_genes.genes.annot
 grep -P "^${gene}\s" ADHD_PGC_2022.genes.annot >> ADHD_PGC_2022_GCST90451108_ADHD_PGC_2022_common_genes.genes.annot
done < GCST90451108_ADHD_PGC_2022_meta_common_genes_ids.txt

# Re-run gene-level analysis for each trait separately
magma --bfile g1000_EUR/g1000_eur --pval GCST90451108_magma_format_sumstats.txt ncol=N --gene-annot GCST90451108_GCST90451108_ADHD_PGC_2022_common_genes.genes.annot --out GCST90451108_GCST90451108_ADHD_PGC_2022
magma --bfile g1000_EUR/g1000_eur --pval ADHD_PGC_2022_magma_format_sumstats.txt ncol=N --gene-annot ADHD_PGC_2022_GCST90451108_ADHD_PGC_2022_common_genes.genes.annot --out ADHD_PGC_2022_GCST90451108_ADHD_PGC_2022

# meta-analyse gene-level analyses
magma --meta raw=GCST90451108_GCST90451108_ADHD_PGC_2022.genes.raw,ADHD_PGC_2022_GCST90451108_ADHD_PGC_2022.genes.raw correlations=GCST90451108_ADHD_PGC_2022_ldsr_intercepts.txt --out GCST90451108_ADHD_PGC_2022_common_genes

# run magma gene-set analysis on meta-analysed results
magma --gene-results GCST90451108_ADHD_PGC_2022_common_genes.genes.raw --set-annot c2.cp.v2024.1.Hs.symbols.gmt.txt --out GCST90451108_ADHD_PGC_2022
