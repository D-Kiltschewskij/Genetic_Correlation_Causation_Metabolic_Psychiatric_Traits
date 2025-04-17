#!/bin/bash
conda activate ldsc
cd </path/to/working/dir>
ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out GCST90451106_ADHD_PGC_2022_ldsc --rg GCST90451106_munged.sumstats.gz,ADHD_PGC_2022_munged.sumstats.gz
