These contain the code for "Dispersion Estimation and Its Effect on Test 
Performance in RNA-seq Data Analysis" by Will Landau and Peng Liu. 

DEPENDENCIES:

1) R version 2.15.3 (2013-03-01). See http://cran.r-project.org/ for details and usage.

2) The following R packages:

  abind version 1.4-0 (hosted on CRAN)
  AMAP.Seq version 1.0 (hosted on CRAN)
  Biobase version 2.18.0 (hosted on Bioconductor)
  clusterGeneration version 1.3.1 (hosted on CRAN)
  DESeq version 1.10.1 (hosted on Bioconductor)
  DSS version 1.0.0 (hosted on Bioconductor)
  edgeR version 3.0.8 (hosted on Bioconductor)
  ggplot2 version 0.9.3.1 (hosted on CRAN)
  hexbin version 1.26.1 (hosted on CRAN)
  iterators version 1.0.6 (hosted on CRAN)
  magic version 1.5-4 (hosted on CRAN)
  MASS version 7.3-23 (hosted on CRAN)
  multicore version 0.1-7 (hosted on CRAN)
  plyr version 1.8 (hosted on CRAN)
  QuasiSeq version 1.0-2 (hosted on CRAN)
  pracma version 1.4.5 (hosted on CRAN)
  reshape2 version 1.2.2 (hosted on CRAN)

The packages hosted on CRAN can be installed by typing the command in R:

  install.packages("packagename")

For help (for example, to find out how to locally install a package), type into R:

  ?install.packages

The packages hosted on Bioconductor can be installed by typing the following commands into R:

  source("http://bioconductor.org/biocLite.R")
  biocLite()
  biocLite("packagename")

For help on installing packages from Bioconductor, type:
  
  source("http://bioconductor.org/biocLite.R")
  ?biocLite


USAGE:

  The commands in example.r demonstrate a basic, step-by-step usage of the simulation code. 
  However, one can run the full code in batch mode in two steps:

1) Configure the makeArgs function in functions.r to set the simulation parameters:

  makeArgs = function(){
    print("Setting simulation parameters.")
    groups = rep(list(c(rep(1,3), rep(2,3)),
                      c(rep(1,3), rep(2,15)),
                      c(rep(1,9), rep(2,9))), times = 2)
    sim.settings = 1:6
  
    lapply(1:length(sim.settings),
           function(cond){list(
             nG = 10000, # number of genes (10000)
             blocksize = 50, # size of each group of correlated pseudo-genes (50)
             nreps = 30, # num. pseudo-datasets per simulation setting (30)
             ncores = 10, # 10
             p0.true = .8,
             sim.setting = sim.settings[cond],
             group = groups[cond][[1]]
           )})
  }

The user has the option to set:

  nG, the number of genes used in each pseudo-dataset.
  blocksize, the size of each group of pseudo-genes with correlated log fold changes
  nreps, the number of pseudo-datasets to generate per simulation setting
  ncores, the number of cores to use. With 10 cores, the workflow takes about a day.
  p0.true, the proportion of equivalently expressed pseudo-genes

2) Modify the path to Rscirpt in run.sh. For example, if Rscirpt is stored in the /usr/local/bin 
directory and run.sh has the line,

  rscript=/apps/R-2.15.3/bin/Rscript 
   
change it to:

  rscript=/usr/local/bin/Rscript

3) In a command line tool like Terminal in Mac OS X, change directories into the directory
containing functions.r. For example:

  cd /myHome/RNASeq_Dispersion_2013/code

4) In the same command line tool, type:

  ./run.sh 

The script takes about a day to run with 10 parallel CPU cores. 
It will create directories "rds" and "fig" to store the following results:

  rds
    args.rds: simulation parameters
    auc.rds: data structure storing the AUC (area under ROC curve) values
    countsI.rds: pseudo-counts for simulation setting I
    countsII.rds: pseudo-counts for simulation setting I
    countsIII.rds: pseudo-counts for simulation setting I
    countsIV.rds: pseudo-counts for simulation setting I
    countsV.rds: pseudo-counts for simulation setting I
    countsVI.rds: pseudo-counts for simulation setting I
    de.rds: array of indicators: 0 means that a pseudo-gene was truly equivalently expressed,
            1 means that a pseudo-gene was truly differentially expressed.
    delta.rds: array of log fold changes, delta. See the simulation section of the main paper.
    phi.rds: array of dispersions
    pvalsI.rds: p-values from the tests for differential expression, simulation setting I
    pvalsII.rds: p-values from the tests for differential expression, simulation setting II
    pvalsIII.rds: p-values from the tests for differential expression, simulation setting II
    pvalsIV.rds: p-values from the tests for differential expression, simulation setting IV
    pvalsV.rds: p-values from the tests for differential expression, simulation setting V
    pvalsVI.rds: p-values from the tests for differential expression, simulation setting VI
    rocI.rds: ROC (receiver operating characteristic) curves for simulation setting I
    rocII.rds: ROC (receiver operating characteristic) curves for simulation setting II
    rocIII.rds: ROC (receiver operating characteristic) curves for simulation setting III
    rocIV.rds: ROC (receiver operating characteristic) curves for simulation setting IV
    rocV.rds: ROC (receiver operating characteristic) curves for simulation setting V
    rocVI.rds: ROC (receiver operating characteristic) curves for simulation setting VI
    sessionInfo.rds: technical info of the current R session: computer architecture, packages, etc.
    true_dispsI.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting I
    true_dispsII.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting II
    true_dispsIII.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting III
    true_dispsIV.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting IV
    true_dispsV.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting V
    true_dispsVI.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting VI
    true_meansI.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting I
    true_meansII.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting II
    true_meansIII.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting III
    true_meansIV.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting IV
    true_meansV.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting V
    true_meansVI.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting VI

  fig
    auc1.pdf: AUC values for simulation setting I
    auc2.pdf: AUC values for simulation setting II
    auc3.pdf: AUC values for simulation setting III
    auc4.pdf: AUC values for simulation setting IV
    auc5.pdf: AUC values for simulation setting V
    auc6.pdf: AUC values for simulation setting VI
    hists.pdf: histograms of means and dispersions of Hammer and Pickrell datasets
               along with similar histograms of example Hammer-generated and Pickrell-generated
               pseudo-datasets
    meandisp_scatter.pdf: mean-dispersion relationships for Hammer, Pickrell, Hammer-generated,
                          and Pickrell-generated data.
    mse.pdf: mean squared error of transformed dispersions
    phi_vs_phiI.pdf: scatterplots of estimated dispersions vs. true dispersions for an example
                     pseudo-dataset from simulation setting I.
    phi_vs_phiII.pdf: scatterplots of estimated dispersions vs. true dispersions for an example
                     pseudo-dataset from simulation setting II.
    phi_vs_phiIII.pdf: scatterplots of estimated dispersions vs. true dispersions for an example
                     pseudo-dataset from simulation setting III.
    phi_vs_phiIV.pdf: scatterplots of estimated dispersions vs. true dispersions for an example
                     pseudo-dataset from simulation setting IV.
    phi_vs_phiV.pdf: scatterplots of estimated dispersions vs. true dispersions for an example
                     pseudo-dataset from simulation setting V.
    phi_vs_phiVI.pdf: scatterplots of estimated dispersions vs. true dispersions for an example
                     pseudo-dataset from simulation setting VI.
    roc.pdf: 2 example ROC curves using tests for differential expression on pseudo-data

  