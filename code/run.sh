#!/bin/bash

# Use the variable, rscript, to store the full path
# of the specific Rscript executable you're using.
# Set rscript=Rscript to use the default version. 
# However, be sure that the version of Rscript you're 
# using is the one corresponding to R version 2.15.3

rscript=/apps/R-2.15.3/bin/Rscript 
#rscript=Rscript

# Actual commands to run the simulations:
nice -19 nohup $rscript -e \
  'source("functions.r"); workflow()' \
  1>> stdout.txt 2>> stderr.txt &