#!/bin/bash
cd </path/to/working/directory>
Rscript CAUSE_wrapper.R --exp GCST90451106_cause_format_sumstats.txt.gz --exp_name GCST90451106 --exp_cols SNP,BETA,SE,A1,A2,P --oc ADHD_PGC_2022_cause_format_sumstats.txt.gz --oc_name ADHD_PGC_2022  --oc_cols SNP,BETA,SE,A1,A2,P --ref g1000/EUR --path . --print_all "No"
