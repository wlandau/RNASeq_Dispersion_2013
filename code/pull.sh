#!/bin/bash

Rscript -e 'if(!file.exists("../rds")) dir.create("../rds")'
scp -P 323 landau@linux11.stat.iastate.edu:"~/RNASeq_Dispersion_2013/rds/*.rds" ../rds
Rscript -e 'source("functions.r"); session(); print.results()'
