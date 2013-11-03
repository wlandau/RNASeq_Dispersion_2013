The files in this directory contain the code for "Dispersion Estimation and Its Effect on Test 
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
  xtable version 1.7-1 (hosted on CRAN)

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


2) In a command line tool like Terminal in Mac OS X, change directories into the directory
containing functions.r. For example:

  cd /myHome/RNASeq_Dispersion_2013/code

3) In the same command line tool, type:

  ./run.sh

Note: the script will create directories "fig" and "rds" and write
the following files:

  fig
    hists_nolegend.pdf: Figure 1 in main paper. Histograms of means and dispersions 
                         of Hammer and Pickrell datasets along with similar histograms 
                         of example Hammer-generated and Pickrell-generated pseudo-datasets
    hists.pdf: same as above except with a legend.
    hists_slides.pdf: same as above except with the default background of ggplot2.

    meandisp_scatter_nolegend.pdf: Figure 2 in main paper. Mean-dispersion relationships 
                                   for Hammer, Pickrell, Hammer-generated, and 
                                   Pickrell-generated data.
    meandisp_scatter.pdf: same as above except with a legend.
    meandisp_scatter_slides.pdf: same as above except with the default background of ggplot2.

    mse.pdf:  Figure 3 in the main paper. Mean squared error of transformed dispersions
    mse_slides.pdf: same as above except with the default background of ggplot2.

    phi_vs_phiI_nolegend.pdf: Scatterplots of estimated dispersions 
                              vs. true dispersions for an example pseudo-dataset from simulation 
                              setting I.
    phi_vs_phiI.pdf: same as above except with a legend.
    phi_vs_phiI_slides.pdf: same as above except with the default background of ggplot2.

    phi_vs_phiII_nolegend.pdf: Figure 4 in the main paper. Scatterplots of estimated dispersions 
                              vs. true dispersions for an example pseudo-dataset from simulation 
                              setting II.
    phi_vs_phiII.pdf: same as above except with a legend.

    phi_vs_phiIII_nolegend.pdf:  Scatterplots of estimated dispersions 
                              vs. true dispersions for an example pseudo-dataset from simulation 
                              setting III.
    phi_vs_phiIII.pdf: same as above except with a legend.
    phi_vs_phiIII_slides.pdf: same as above except with the default background of ggplot2.

    phi_vs_phiIV_nolegend.pdf:  Scatterplots of estimated dispersions 
                              vs. true dispersions for an example pseudo-dataset from simulation 
                              setting IV.
    phi_vs_phiIV.pdf: same as above except with a legend.
    phi_vs_phiIV_slides.pdf: same as above except with the default background of ggplot2.

    phi_vs_phiV_nolegend.pdf: Figure 5 in the main paper. Scatterplots of estimated dispersions 
                              vs. true dispersions for an example pseudo-dataset from simulation 
                              setting V.
    phi_vs_phiV.pdf: same as above except with a legend.
    phi_vs_phiV_slides.pdf: same as above except with the default background of ggplot2.

    phi_vs_phiVI_nolegend.pdf:  Scatterplots of estimated dispersions 
                              vs. true dispersions for an example pseudo-dataset from simulation 
                              setting VI.
    phi_vs_phiVI.pdf: same as above except with a legend.
    phi_vs_phiVI_slides.pdf: same as above except with the default background of ggplot2.

    auc1.pdf: Figure 6 in main paper. AUC (area under ROC curve) values for simulation setting I
    auc1_slides.pdf: same as above except with the default background of ggplot2.

    auc2.pdf: Figure 7 in main paper. AUC values for simulation setting II
    auc2_slides.pdf: same as above except with the default background of ggplot2.

    auc3.pdf: Figure 8 in main paper. AUC values for simulation setting III
    auc3_slides.pdf: same as above except with the default background of ggplot2.

    auc4.pdf: Figure 9 in the main paper. AUC values for simulation setting IV
    auc4_slides.pdf: same as above except with the default background of ggplot2.

    auc5.pdf: Figure 10 in the main paper. AUC values for simulation setting V
    auc5_slides.pdf: same as above except with the default background of ggplot2.

    auc6.pdf: Figure 11 in the main paper. AUC values for simulation setting VI
    auc6_slides.pdf: same as above except with the default background of ggplot2.

    roc_nolegend.pdf: 2 example ROC curves using tests for differential expression on pseudo-data.
    roc.pdf: same as above except with a legend.
    roc_slides.pdf: same as above except with the default background of ggplot2.

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

  

  