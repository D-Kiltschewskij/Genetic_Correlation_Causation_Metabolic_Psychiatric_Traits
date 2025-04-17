#!/bin/bash
cd </path/to/working/dir>
Rscript LCV_wrapper.R --ldscores all.chr.l2.ldscore --path . --phenotype_1 GCST90451106 --sumstats_1 GCST90451106_munged.sumstats.gz --phenotype_2 ADHD_PGC_2022 --sumstats_2 ADHD_PGC_2022_munged.sumstats.gz
