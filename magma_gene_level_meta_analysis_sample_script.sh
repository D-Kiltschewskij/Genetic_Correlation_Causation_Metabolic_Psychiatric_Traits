#!/bin/bash

conda activate ldsc

cd </path/to/working/directory>

# Estimate genetic correlation
ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out GCST90451108_GCST90451108_ldsc --rg GCST90451108_munged.sumstats.gz,GCST90451108_munged.sumstats.gz
ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out GCST90451108_ADHD_PGC_2022_ldsc --rg GCST90451108_munged.sumstats.gz,ADHD_PGC_2022_munged.sumstats.gz
ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ADHD_PGC_2022_GCST90451108_ldsc --rg ADHD_PGC_2022_munged.sumstats.gz,GCST90451108_munged.sumstats.gz
ldsc.py --ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ --out ADHD_PGC_2022_ADHD_PGC_2022_ldsc --rg ADHD_PGC_2022_munged.sumstats.gz,/ADHD_PGC_2022_munged.sumstats.gz

# Assemble matrix of ldsc genetic covariance intercepts
EXPEXP=$(grep "Intercept" GCST90451108_GCST90451108_ldsc.log | tail -n 1 | awk -F " " '{print $2}')
EXPOUT=$(grep "Intercept" GCST90451108_ADHD_PGC_2022_ldsc.log | tail -n 1 | awk -F " " '{print $2}')
OUTEXP=$(grep "Intercept" ADHD_PGC_2022_GCST90451108_ldsc.log | tail -n 1 | awk -F " " '{print $2}')
OUTOUT=$(grep "Intercept" ADHD_PGC_2022_ADHD_PGC_2022_ldsc.log | tail -n 1 | awk -F " " '{print $2}')
echo "${EXPEXP} ${OUTEXP}" > GCST90451108_ADHD_PGC_2022_ldsr_intercepts.txt
echo "${EXPOUT} ${OUTOUT}" >> GCST90451108_ADHD_PGC_2022_ldsr_intercepts.txt

# Run MAGMA meta-analysis
magma --meta genes=GCST90451108.genes.out,ADHD_PGC_2022.genes.out correlations=GCST90451108_ADHD_PGC_2022_ldsr_intercepts.txt --out GCST90451108_ADHD_PGC_2022_meta
head -1 GCST90451108_ADHD_PGC_2022_meta.genes.out > GCST90451108_ADHD_PGC_2022_meta_common_genes.txt
awk -F " " '$8 > 1' GCST90451108_ADHD_PGC_2022_meta.genes.out > GCST90451108_ADHD_PGC_2022_meta_common_genes.txt
