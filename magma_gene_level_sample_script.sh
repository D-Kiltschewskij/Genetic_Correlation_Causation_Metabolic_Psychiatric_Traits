#!/bin/bash
cd </path/to/working/directory>

magma --annotate window=5,1.5 --snp-loc ADHD_PGC_2022_magma_format_sumstats.txt --gene-loc NCBI37_MHCRMV_clean.3.gene.loc --out ADHD_PGC_2022

magma --bfile g1000_EUR/g1000_eur --pval ADHD_PGC_2022_magma_format_sumstats.txt ncol=N --gene-annot ADHD_PGC_2022.genes.annot --out ADHD_PGC_2022
