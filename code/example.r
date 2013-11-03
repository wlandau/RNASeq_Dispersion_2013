# You must be in the directory containing functions.r.
# For example:
setwd("~/RNASeq_Dispersion_2013/code")

# Load everything in functions.r
source("functions.r")

# Use the initialize() function to load all
# required libraries and save session info.
initialize()

# Now, I define some global parameters 
# in a list to pass into the utility functions.
arg = list(
  nG = 100, # number of genes in each pseudo-dataset
  blocksize = 2, # number of genes in each block of correlated genes
  nreps = 1, # number of pseudo-datasets to generate
  p0.true = .2, # percent of truly null genes in the pseudo-datasets
  ncores = 1, # number of CPU cores to use in the calculations 
  sim.setting = 1, # simulation setting
  tst.setting = 1, # test setting (possibly unnecessary)
  dsp.setting = "A", # dispersion estimation setting (possibly unnecessary)
  group = c(1,1,1,2,2,2) # assignment of libraries to treatment groups
)

# compute input to simulations using
# data/montpick_eset.RData and
# data/hammer_eset.RData
computeInput()

# simulate a count dataset using information saved by computeInput()
sim = sim.counts(arg)

# sim is now a list that contains the pseudo-dataset,
# a record of which genes are truly DE, and true 
# dispersions. sim will eventually contain results
# from analyses on the pseudo-dataset.
str(sim)

# I can use edgeR to estimate tagwise dispersions on the 
# count dataset
sim = wqCML.dispersions(sim)
str(sim)
head(sim$phi)

# I can apply the exact test in edgeR to the data using
# both the true dispersions and the edgeR qCML dispersions
sim = edgeR.test(sim)
str(sim)
head(sim$pvals)