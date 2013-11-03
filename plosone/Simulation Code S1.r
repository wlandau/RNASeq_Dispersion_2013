# Author: Will Landau
# Organization: Iowa State University
# Email: landau@iastate.edu

## README ##

		# This file contains the simulation code for "Dispersion Estimation and Its Effect on Test 
		# Performance in RNA-seq Data Analysis" by Will Landau and Peng Liu. 
		
		# DEPENDENCIES:
		
		# 1) R version 2.15.3 (2013-03-01). See http://cran.r-project.org/ for details and usage.
		
		# 2) The following R packages:
		
		  # abind version 1.4-0 (hosted on CRAN)
		  # AMAP.Seq version 1.0 (hosted on CRAN)
		  # Biobase version 2.18.0 (hosted on Bioconductor)
		  # clusterGeneration version 1.3.1 (hosted on CRAN)
		  # DESeq version 1.10.1 (hosted on Bioconductor)
		  # DSS version 1.0.0 (hosted on Bioconductor)
		  # edgeR version 3.0.8 (hosted on Bioconductor)
		  # ggplot2 version 0.9.3.1 (hosted on CRAN)
		  # hexbin version 1.26.1 (hosted on CRAN)
		  # iterators version 1.0.6 (hosted on CRAN)
		  # magic version 1.5-4 (hosted on CRAN)
		  # MASS version 7.3-23 (hosted on CRAN)
		  # multicore version 0.1-7 (hosted on CRAN)
		  # plyr version 1.8 (hosted on CRAN)
		  # QuasiSeq version 1.0-2 (hosted on CRAN)
		  # pracma version 1.4.5 (hosted on CRAN)
		  # reshape2 version 1.2.2 (hosted on CRAN)
		  # xtable version 1.7-1 (hosted on CRAN)
		
		# The packages hosted on CRAN can be installed by typing the command in R:
		
		  # install.packages("packagename")
		
		# For help (for example, to find out how to locally install a package), type into R:
		
		  # ?install.packages
		
		# The packages hosted on Bioconductor can be installed by typing the following commands into R:
		
		  # source("http://bioconductor.org/biocLite.R")
		  # biocLite()
		  # biocLite("packagename")
		
		# For help on installing packages from Bioconductor, type:
		  
		  # source("http://bioconductor.org/biocLite.R")
		  # ?biocLite
		
		
		# USAGE:
		
		  # The commands in example.r demonstrate a basic, step-by-step usage of the simulation code. 
		  # However, one can run the full code in batch mode in two steps:
		
		# 1) Configure the makeArgs function in functions.r to set the simulation parameters:
		
		  # makeArgs = function(){
		    # print("Setting simulation parameters.")
		    # groups = rep(list(c(rep(1,3), rep(2,3)),
		                      # c(rep(1,3), rep(2,15)),
		                      # c(rep(1,9), rep(2,9))), times = 2)
		    # sim.settings = 1:6
		  
		    # lapply(1:length(sim.settings),
		           # function(cond){list(
		             # nG = 10000, # number of genes (10000)
		             # blocksize = 50, # size of each group of correlated pseudo-genes (50)
		             # nreps = 30, # num. pseudo-datasets per simulation setting (30)
		             # ncores = 10, # 10
		             # p0.true = .8,
		             # sim.setting = sim.settings[cond],
		             # group = groups[cond][[1]]
		           # )})
		  # }
		
		# The user has the option to set:
		
		  # nG, the number of genes used in each pseudo-dataset.
		  # blocksize, the size of each group of pseudo-genes with correlated log fold changes
		  # nreps, the number of pseudo-datasets to generate per simulation setting
		  # ncores, the number of cores to use. With 10 cores, the workflow takes about a day.
		  # p0.true, the proportion of equivalently expressed pseudo-genes
		
		
		# 2) In a command line tool like Terminal in Mac OS X, change directories into the directory
		# containing functions.r. For example:
		
		  # cd /myHome/RNASeq_Dispersion_2013/code
		
		# 3) In the same command line tool, type:
		
		  # ./run.sh
		
		# Note: the script will create directories "fig" and "rds" and write
		# the following files:
		
		  # fig
		    # hists_nolegend.pdf: Figure 1 in main paper. Histograms of means and dispersions 
		                         # of Hammer and Pickrell datasets along with similar histograms 
		                         # of example Hammer-generated and Pickrell-generated pseudo-datasets
		    # hists.pdf: same as above except with a legend.
		    # hists_slides.pdf: same as above except with the default background of ggplot2.
		
		    # meandisp_scatter_nolegend.pdf: Figure 2 in main paper. Mean-dispersion relationships 
		                                   # for Hammer, Pickrell, Hammer-generated, and 
		                                   # Pickrell-generated data.
		    # meandisp_scatter.pdf: same as above except with a legend.
		    # meandisp_scatter_slides.pdf: same as above except with the default background of ggplot2.
		
		    # mse.pdf:  Figure 3 in the main paper. Mean squared error of transformed dispersions
		    # mse_slides.pdf: same as above except with the default background of ggplot2.
		
		    # phi_vs_phiI_nolegend.pdf: Scatterplots of estimated dispersions 
		                              # vs. true dispersions for an example pseudo-dataset from simulation 
		                              # setting I.
		    # phi_vs_phiI.pdf: same as above except with a legend.
		    # phi_vs_phiI_slides.pdf: same as above except with the default background of ggplot2.
		
		    # phi_vs_phiII_nolegend.pdf: Figure 4 in the main paper. Scatterplots of estimated dispersions 
		                              # vs. true dispersions for an example pseudo-dataset from simulation 
		                              # setting II.
		    # phi_vs_phiII.pdf: same as above except with a legend.
		
		    # phi_vs_phiIII_nolegend.pdf:  Scatterplots of estimated dispersions 
		                              # vs. true dispersions for an example pseudo-dataset from simulation 
		                              # setting III.
		    # phi_vs_phiIII.pdf: same as above except with a legend.
		    # phi_vs_phiIII_slides.pdf: same as above except with the default background of ggplot2.
		
		    # phi_vs_phiIV_nolegend.pdf:  Scatterplots of estimated dispersions 
		                              # vs. true dispersions for an example pseudo-dataset from simulation 
		                              # setting IV.
		    # phi_vs_phiIV.pdf: same as above except with a legend.
		    # phi_vs_phiIV_slides.pdf: same as above except with the default background of ggplot2.
		
		    # phi_vs_phiV_nolegend.pdf: Figure 5 in the main paper. Scatterplots of estimated dispersions 
		                              # vs. true dispersions for an example pseudo-dataset from simulation 
		                              # setting V.
		    # phi_vs_phiV.pdf: same as above except with a legend.
		    # phi_vs_phiV_slides.pdf: same as above except with the default background of ggplot2.
		
		    # phi_vs_phiVI_nolegend.pdf:  Scatterplots of estimated dispersions 
		                              # vs. true dispersions for an example pseudo-dataset from simulation 
		                              # setting VI.
		    # phi_vs_phiVI.pdf: same as above except with a legend.
		    # phi_vs_phiVI_slides.pdf: same as above except with the default background of ggplot2.
		
		    # auc1.pdf: Figure 6 in main paper. AUC (area under ROC curve) values for simulation setting I
		    # auc1_slides.pdf: same as above except with the default background of ggplot2.
		
		    # auc2.pdf: Figure 7 in main paper. AUC values for simulation setting II
		    # auc2_slides.pdf: same as above except with the default background of ggplot2.
		
		    # auc3.pdf: Figure 8 in main paper. AUC values for simulation setting III
		    # auc3_slides.pdf: same as above except with the default background of ggplot2.
		
		    # auc4.pdf: Figure 9 in the main paper. AUC values for simulation setting IV
		    # auc4_slides.pdf: same as above except with the default background of ggplot2.
		
		    # auc5.pdf: Figure 10 in the main paper. AUC values for simulation setting V
		    # auc5_slides.pdf: same as above except with the default background of ggplot2.
		
		    # auc6.pdf: Figure 11 in the main paper. AUC values for simulation setting VI
		    # auc6_slides.pdf: same as above except with the default background of ggplot2.
		
		    # roc_nolegend.pdf: 2 example ROC curves using tests for differential expression on pseudo-data.
		    # roc.pdf: same as above except with a legend.
		    # roc_slides.pdf: same as above except with the default background of ggplot2.
		
		  # rds
		    # args.rds: simulation parameters
		    # auc.rds: data structure storing the AUC (area under ROC curve) values
		    # countsI.rds: pseudo-counts for simulation setting I
		    # countsII.rds: pseudo-counts for simulation setting I
		    # countsIII.rds: pseudo-counts for simulation setting I
		    # countsIV.rds: pseudo-counts for simulation setting I
		    # countsV.rds: pseudo-counts for simulation setting I
		    # countsVI.rds: pseudo-counts for simulation setting I
		    # de.rds: array of indicators: 0 means that a pseudo-gene was truly equivalently expressed,
		            # 1 means that a pseudo-gene was truly differentially expressed.
		    # delta.rds: array of log fold changes, delta. See the simulation section of the main paper.
		    # phi.rds: array of dispersions
		    # pvalsI.rds: p-values from the tests for differential expression, simulation setting I
		    # pvalsII.rds: p-values from the tests for differential expression, simulation setting II
		    # pvalsIII.rds: p-values from the tests for differential expression, simulation setting II
		    # pvalsIV.rds: p-values from the tests for differential expression, simulation setting IV
		    # pvalsV.rds: p-values from the tests for differential expression, simulation setting V
		    # pvalsVI.rds: p-values from the tests for differential expression, simulation setting VI
		    # rocI.rds: ROC (receiver operating characteristic) curves for simulation setting I
		    # rocII.rds: ROC (receiver operating characteristic) curves for simulation setting II
		    # rocIII.rds: ROC (receiver operating characteristic) curves for simulation setting III
		    # rocIV.rds: ROC (receiver operating characteristic) curves for simulation setting IV
		    # rocV.rds: ROC (receiver operating characteristic) curves for simulation setting V
		    # rocVI.rds: ROC (receiver operating characteristic) curves for simulation setting VI
		    # sessionInfo.rds: technical info of the current R session: computer architecture, packages, etc.
		    # true_dispsI.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting I
		    # true_dispsII.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting II
		    # true_dispsIII.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting III
		    # true_dispsIV.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting IV
		    # true_dispsV.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting V
		    # true_dispsVI.rds: dispersions phihat_g used to simulate pseudo-data for simulation setting VI
		    # true_meansI.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting I
		    # true_meansII.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting II
		    # true_meansIII.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting III
		    # true_meansIV.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting IV
		    # true_meansV.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting V
		    # true_meansVI.rds: geometric means ybar_g. used to simulate pseudo-data for simulation setting VI
		
		
		
		# ## EXAMPLE USAGE ##
		# # You must be in the directory containing functions.r.
		# # For example:
		# setwd("~/RNASeq_Dispersion_2013/code")
		
		# # Load everything in functions.r
		# source("functions.r")
		
		# # Use the initialize() function to load all
		# # required libraries and save session info.
		# initialize()
		
		# # Now, I define some global parameters 
		# # in a list to pass into the utility functions.
		# arg = list(
		  # nG = 100, # number of genes in each pseudo-dataset
		  # blocksize = 2, # number of genes in each block of correlated genes
		  # nreps = 1, # number of pseudo-datasets to generate
		  # p0.true = .2, # percent of truly null genes in the pseudo-datasets
		  # ncores = 1, # number of CPU cores to use in the calculations 
		  # sim.setting = 1, # simulation setting
		  # tst.setting = 1, # test setting (possibly unnecessary)
		  # dsp.setting = "A", # dispersion estimation setting (possibly unnecessary)
		  # group = c(1,1,1,2,2,2) # assignment of libraries to treatment groups
		# )
		
		# # compute input to simulations using
		# # data/montpick_eset.RData and
		# # data/hammer_eset.RData
		# computeInput()
		
		# # simulate a count dataset using information saved by computeInput()
		# sim = sim.counts(arg)
		
		# # sim is now a list that contains the pseudo-dataset,
		# # a record of which genes are truly DE, and true 
		# # dispersions. sim will eventually contain results
		# # from analyses on the pseudo-dataset.
		# str(sim)
		
		# # I can use edgeR to estimate tagwise dispersions on the 
		# # count dataset
		# sim = wqCML.dispersions(sim)
		# str(sim)
		# head(sim$phi)
		
		# # I can apply the exact test in edgeR to the data using
		# # both the true dispersions and the edgeR qCML dispersions
		# sim = edgeR.test(sim)
		# str(sim)
		# head(sim$pvals)
		
		
		# ## SHELL SCRIPT FOR GENERATING THE FIGURES IN THE PAPER ##
		
		# #!/bin/bash
		# #Generate figures for the manuscript using the program's output in ../fig.
		
		# outdir="../plosone"
		# indir="../fig"
		
		# if [[ ! -d  $outdir ]];
		# then
		   # mkdir $outdir
		# fi
		
		# # Convert pdfs to tifs ImageMagick
		
		# in=(hists_nolegend.pdf meandisp_scatter_nolegend.pdf mse.pdf phi_vs_phi_II_nolegend.pdf phi_vs_phi_V_nolegend.pdf auc1.pdf auc2.pdf auc3.pdf auc4.pdf auc5.pdf auc6.pdf)
		
		# for (( i=0; i<${#in[@]} ; i++ ))
		# do
		  # convert -strip -units PixelsPerInch -density 300 -resample 300 -alpha off -colorspace RGB -depth 8 -trim -bordercolor white -border 1% -resize '2049x2758>' -resize '980x980<' +repage -compress lzw $indir/${in[$i]} $outdir/Figure$[$i+1].tiff
		# done


## FUNCTIONS ##

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

err = function(expr, sim = NULL, info = NULL){
  tryCatch(expr, error = function(e){
    write(paste(e), stderr())
    write(paste(traceback()), stderr())
    filename = paste("../rds/badsim", 
                     substr(as.character(round(runif(1), 6)), 3, 9), 
                     ".rds", sep = "")
    sim$info = info
    saveRDS(sim, filename)
    stop()
  })  
}

sup = function(x){
  suppressMessages(suppressWarnings(x))
}

session = function(){
  print("Loading libraries.")
  
  sup(library(abind))
  sup(library(AMAP.Seq))
  sup(library(clusterGeneration))
  sup(library(DESeq))
  sup(library(DSS))
  sup(library(edgeR))
  sup(library(ggplot2))
  sup(library(hexbin))
  sup(library(iterators))
  sup(library(magic))
  sup(library(MASS))
  sup(library(multicore))
  sup(library(plyr))
  sup(library(QuasiSeq))
  sup(library(pracma))
  sup(library(reshape2))
}

initialize = function(){
  session()
    
  if(!file.exists("progress.txt")) 
    file.create("progress.txt")
  
  if(!file.exists("../rds"))
    dir.create("../rds")
  
  if(!file.exists("../fig"))
    dir.create("../fig")
  
  seed = 10
  set.seed(seed)
  args = makeArgs()
  
  session.info=list(sessionInfo=sessionInfo(), 
                    Sys.info=Sys.info(),
                    Sys.getenv=Sys.getenv(names=TRUE), 
                    RNGkind=RNGkind(), 
                    .libPaths=.libPaths(), 
                    date=date(), 
                    wd=getwd(),
                    packInfo = installed.packages(), 
                    simulation.parameters = list(seed = seed, 
                                                 args = args))
  
  saveRDS(args, "../rds/args.rds")
  saveRDS(session.info, "../rds/sessionInfo.rds")
}

geom.mean = function(row){
  row[row == 0] = 0.1
  if(length(row) != 0){
    return(exp(sum(log(row))/length(row)))
  } else{
    return(0)
  }
}

geom.rowMeans = function(m){
  apply(m, 1, geom.mean)
}

pickrell = function(){
  load(url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData"))
#  load("../data/montpick_eset.RData")
  counts = exprs(montpick.eset)
  counts = counts[, phenoData(montpick.eset)$study == "Pickrell"]
  counts = counts[rowSums(counts) > 0, ]
  group = as.factor(rep(1, dim(counts)[2]))
  list(counts = counts, group = group)
}

hammer = function(){
  load(url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/hammer_eset.RData"))
#  load("../data/hammer_eset.RData")
  counts = exprs(hammer.eset)
  counts = counts[rowSums(counts) > 0,]
  
  group = as.vector(phenoData(hammer.eset)$protocol)
  group[group == "control"] = "1"
  group[group == "L5 SNL"] = "2"
  group = as.factor(as.numeric(group))
  
  list(counts = counts, group = group)
}

computeInput = function(){
  print("Computing input.")
  sim = pickrell()

  d <- newCountDataSet(sim$counts, sim$group)
  d <- estimateSizeFactors(d)
  norm = pData(d)$sizeFactor

  dsp = dispersion.nb.QL(sim$counts, 
                         norm, 
                         sim$group)$dispersion

  rmeans = geom.rowMeans(sim$counts)
  
  for(setting in 1:6){
    if(setting == 4){
      sim = hammer()

      d <- newCountDataSet(sim$counts, sim$group)
      d <- estimateSizeFactors(d)
      norm = pData(d)$sizeFactor

      dsp = dispersion.nb.QL(sim$counts, 
                             norm, 
                             sim$group)$dispersion

      rmeans = geom.rowMeans(sim$counts)
    }
      
    saveRDS(rmeans,
            paste("../rds/true_means", as.roman(setting), ".rds", sep=""))
    
    saveRDS(dsp,
            paste("../rds/true_disps", as.roman(setting), ".rds", sep=""))
  }
}

next.g = function(env){
  tryCatch(nextElem(env$gene.iter), error=function(e){
    print("WARNING: recycling genes.")
    write("WARNING: recycling genes.", stderr())
    env$gene.iter = iter(sample(1:env$gene.iter$length))
    return(nextElem(env$gene.iter))
  })
}

next.genes = function(n, env){
  replicate(n, next.g(env))
}

makeSigma = function(arg){  
  n = round(arg$nG - arg$nG * arg$p0.true)
  if(arg$blocksize == 1)
    return(diag(rep(1, n)))
  
  if(!is.finite(arg$blocksize) || is.null(arg$blocksize) || arg$blocksize < 1)
    arg$blocksize = 1
  
  arg$blocksize = round(arg$blocksize)
  dims = rep(arg$blocksize, n %/% arg$blocksize)
  if(s <- n - sum(dims))
    dims = c(dims, s)
  
  mats = lapply(as.list(dims), rcorrmatrix)
  do.call("adiag", mats)
}

resample.delta = function(delta, de, h, Sigma){
  delta.de = delta[de]
  resample = h[de]
  
  if (all(!resample)){
    return(delta)
  } else if(all(resample)){
    conditional.mean = rep(0, length(resample))
    conditional.variance = Sigma
  } else {
    s11 = matrix(Sigma[resample, resample], nrow = sum(resample), ncol = sum(resample))
    s12 = matrix(Sigma[resample, !resample], nrow = sum(resample), ncol = sum(!resample))
    s21 = matrix(Sigma[!resample, resample], nrow = sum(!resample), ncol = sum(resample))
    s22 = matrix(Sigma[!resample, !resample], nrow = sum(!resample), ncol = sum(!resample))
    s22inv = solve(s22)
    
    conditional.mean = s12 %*% s22inv %*% delta.de[!resample]
    conditional.variance = s11 - s12 %*% s22inv %*% s21
  } 
  
  delta.de[resample] = mvrnorm(1, conditional.mean, conditional.variance)
  delta[de] = delta.de
  delta
}

sim.counts = function(arg){
  print("Simulating a pseudo-dataset.")
  
  true_means = readRDS(paste("../rds/true_means", 
                             as.roman(arg$sim.setting), ".rds", sep=""))
  true_disps = readRDS(paste("../rds/true_disps", 
                             as.roman(arg$sim.setting), ".rds", sep=""))
  
  env = new.env()
  assign("gene.iter",iter(sample(1:length(true_means))), envir=env)
  
  de = c(rep(F, round(arg$nG * arg$p0.true)), 
         rep(T, round(arg$nG - arg$nG * arg$p0.true)))
  de = de[sample.int(length(de))] 
  
  # h = vector indicating which pseudo-genes to be re-simulated
  h = rep(T, arg$nG) 
  counts <- phi <- matrix(0, nrow = arg$nG, ncol = length(arg$group))
  
  # delta = log fold changes. delta[de] ~ multivariate normal (0, Sigma)
  delta <- selected_genes <- selected_means <- rep(0, arg$nG)
  lambda = matrix(0, nrow = arg$nG, ncol = length(arg$group))
  prms = array(0, dim = c(arg$nG, length(arg$group), 2))
  
  print("  Creating the covariance matrix of the pseudo-data.")
  Sigma = makeSigma(arg)
  
  print("  Drawing pseudo-counts.")
  while(any(h)){
    delta = resample.delta(delta, de, h, Sigma)
    
    selected_genes[h] = next.genes(sum(h), env)
    selected_means[h] = true_means[selected_genes[h]]
    
    lambda[h,] = outer(delta[h], arg$group, function(x, y){
      exp((-1)^y * (1/2) * x)
    })
    lambda[h,] = sweep(matrix(lambda[h,], nrow = sum(h), ncol = length(arg$group)), 
                       1, selected_means[h], "*")
    
    phi[h,] = replicate(length(arg$group), true_disps[selected_genes[h]]) 
    prms[h,,] = abind(matrix(lambda[h,], ncol = length(arg$group)), 
                 matrix(phi[h,], ncol = length(arg$group)), 
                 along = 3)
    
    counts[h,] = apply(array(prms[h,,], dim = c(sum(h), length(arg$group), 2)), 
                       c(1,2), 
                       function(x){rnegbin(1, x[1], 1/x[2])})
    h = (rowSums(counts) == 0)
  }
   
  if(any(rowSums(counts) <= 0))
    err(1 <- 0, 
        sim = list(arg = arg, counts = counts, de = de, h = h,
                          delta = delta, Sigma = Sigma),
        info = "Failed to simulate data: some genes are not expressed.")
  
  nms = paste("gene", 1:arg$nG, sep="")
  names(delta) = nms
  rownames(counts) = nms
  colnames(counts) = paste("lib", 1:ncol(counts), sep="")
  
  phi = matrix(phi[,1], ncol = 1)
  rownames(phi) = nms
  colnames(phi) = "True"
  
  names(de) = nms

  d <- newCountDataSet(counts, arg$group)
  d <- estimateSizeFactors(d)
  norm = pData(d)$sizeFactor  # Normalization factors computed with the
                              # method by Anders and Huber (2010).
  
  gm = prod(norm)^(1/length(norm))
  md = median(apply(counts, 2, sum))

  size = floor(norm * md /gm) # Adjusted library sizes for edgeR functions. Setting
                              # norm.factors = size in the DGEList() function adapts
                              # Anders and Huber normalization method to the 
                              # differential expression analysis methods in edgeR.

  if(length(unique(size)) == 1)
    size[1] = size[1] + 1  

  list(counts = counts, 
       group = arg$group, 
       phi = phi, 
       de = de, 
       delta = delta,
       size = size,
       norm = norm, 
       nG = arg$nG,
       sim.setting = arg$sim.setting)
}

QL.dispersions = function(sim, dsp.setting){
  print(paste("  ", dsp.setting))
  nms = colnames(sim$phi)
  sim$phi = cbind(sim$phi, 
                  dispersion.nb.QL(sim$counts, 
                                   sim$norm, 
                                   sim$group)$dispersion)
  colnames(sim$phi) = c(nms, dsp.setting)
  return(sim)
}

DSS.dispersions = function(sim, 
                           dsp.setting = "", 
                           trend = F){
  print(paste("  ", dsp.setting))
  cts = sim$counts
  colnames(cts) = NULL

  d = newSeqCountSet(cts, sim$group)
  attr(d, "normalizationFactor") = sim$norm
  d = estDispersion(d, trend)
  
  nms = colnames(sim$phi)
  sim$phi = cbind(sim$phi, dispersion(d))
  colnames(sim$phi) = c(nms, dsp.setting)
  return(sim)
}

wqCML.dispersions = function(sim, common = F, dsp.setting = ""){
  print(paste("  ", 
              dsp.setting))

  d = DGEList(group = sim$group,
              lib.size = sim$size,
              counts = sim$counts)
  d = estimateCommonDisp(d)
  
  if(!common)
    d = estimateTagwiseDisp(d)
  
  nms = colnames(sim$phi)
  
  if(!common)
    sim$phi = cbind(sim$phi, d$tagwise.dispersion)
  else
    sim$phi = cbind(sim$phi, rep(d$common.dispersion, dim(sim$phi)[1]))
  
  colnames(sim$phi) = c(nms, dsp.setting)
  return(sim)
}

APL.dispersions = function(sim, 
                           dsp.setting = "", 
                           method){
  print(paste("  ", dsp.setting))

  d = DGEList(group = sim$group,
              lib.size = sim$size,
              counts = sim$counts)
  design = model.matrix(~ sim$group)
  d = estimateGLMCommonDisp(y = d, design = design)
  
  if(method == "Common"){
    phi.new = rep(d$common.dispersion, sim$nG)
  } else if(method == "Trended"){
    d = estimateGLMTrendedDisp(y = d, design = design)
    phi.new = d$trended.dispersion
  } else if(method == "Tagwise"){
    d = estimateGLMTagwiseDisp(y = d, design = design)
    phi.new = d$tagwise.dispersion
  }
  
  nms = colnames(sim$phi)
  sim$phi = cbind(sim$phi, phi.new)
  colnames(sim$phi) = c(nms, dsp.setting)
  return(sim)
}

DESeq.dispersions = function(sim, dsp.setting = "", 
                             sharingMode = "maximum", 
                             fitType = "local"){
  print(paste("  ", dsp.setting))
  d <- newCountDataSet(sim$counts, sim$group)

#  d <- estimateSizeFactors(d)
  pData(d)$sizeFactor = sim$norm

  d <- suppressWarnings(estimateDispersions(d, method = "pooled",
                                            fitType = fitType,
                                            sharingMode = sharingMode))
  
  nms = colnames(sim$phi)
  sim$phi = cbind(sim$phi, fData(d)$disp_pooled)
  colnames(sim$phi) = c(nms, dsp.setting)
  return(sim)
}

results2roc = function(results){

  results = results[order(results$pvals),]
  results$pvals[!is.finite(results$pvals)] = 1

  pvals.de <- results$pvals[results$de == 1]
  pvals.ee <- results$pvals[results$de == 0]
  cutoffs <- 1:1000 * 1e-3
  fpr <- tpr <- rep(1, 1000)
  
  for(i in 1:1000){
    fpr[i] = sum(pvals.ee <= cutoffs[i])/length(pvals.ee)
    tpr[i] = sum(pvals.de <= cutoffs[i])/length(pvals.de) 
  }
  
  stepf <- stepfun(x = fpr, y = c(0, tpr))
  fpr <- 0:999 / 999
  tpr <- stepf(fpr)
  
  list(tpr = tpr, 
       auc = trapz(x = fpr[fpr < 0.1], 
                   y = tpr[fpr < 0.1]))
}

update.sim = function(sim, tst.setting, 
                      roc, auc, pvals, 
                      pvals.de, pvals.ee){
  nms = c(dimnames(sim$roc)[[3]], 
          paste("test", tst.setting, sep = ""))
  
  colnames(roc) = colnames(sim$phi)
  sim$roc = abind(sim$roc, roc, along = 3)
  dimnames(sim$roc)[[3]] = nms
  
  sim$auc = rbind(sim$auc, auc)
  colnames(sim$auc) = colnames(sim$phi)
  rownames(sim$auc) = nms
  
  colnames(pvals) = colnames(sim$phi)
  sim$pvals = abind(sim$pvals, pvals, along = 3)
  dimnames(sim$pvals)[[1]] = rownames(sim$counts)
  dimnames(sim$pvals)[[3]] = nms
  
  colnames(pvals.de) = colnames(sim$phi)
  sim$pvals.de = abind(sim$pvals.de, pvals.de, along = 3)
  dimnames(sim$pvals.de)[[1]] = rownames(sim$counts)[sim$de == 1]
  dimnames(sim$pvals.de)[[3]] = nms
  
  colnames(pvals.ee) = colnames(sim$phi)
  sim$pvals.ee = abind(sim$pvals.ee, pvals.ee, along = 3)
  dimnames(sim$pvals.ee)[[1]] = rownames(sim$counts)[sim$de == 0]
  dimnames(sim$pvals.ee)[[3]] = nms
  
  sim
}

edgeR.test = function(sim){
  print("Doing edgeR tests.")
  roc <- auc <- pvals <- pvals.de <- pvals.ee <- NULL
  
  for(i in 1:dim(sim$phi)[2]){
    print(paste("  Testing with dispersions", colnames(sim$phi)[i]))
    d = DGEList(group = sim$group, 
                lib.size = sim$size,
                counts = sim$counts)

    d$tagwise.dispersion <- sim$phi[,i]
    
    results = exactTest(d, dispersion = "tagwise")$table
    results$de <- sim$de
    results$pvals = results$PValue
    
    roc.raw = results2roc(results)
    roc = cbind(roc, roc.raw$tpr)
    auc = c(auc, roc.raw$auc)
    
    pvals = cbind(pvals, results$pvals)
    pvals.de = cbind(pvals.de, results$pvals[sim$de == 1])
    pvals.ee = cbind(pvals.ee, results$pvals[sim$de == 0])
  }
  
  update.sim(sim, tst.setting = 1, roc, auc, 
             pvals, pvals.de, pvals.ee)
}

DESeq.test = function(sim){
  print("Doing DESeq tests.")
  roc <- auc <- pvals <- pvals.de <- pvals.ee <- NULL
  
  for(i in 1:dim(sim$phi)[2]){
    print(paste("  Testing with dispersions", colnames(sim$phi)[i]))
    d = newCountDataSet(sim$counts, sim$group)

 #   d <- estimateSizeFactors(d)
    pData(d)$sizeFactor = sim$norm

    d = estimateDispersions(d, fitType = "local")
    fData(d) = data.frame(disp_pooled = sim$phi[,i])
    
    results <- nbinomTest( d, "1", "2")
    results$de = sim$de
    results$pvals = results$padj
    
    roc.raw = results2roc(results)
    roc = cbind(roc, roc.raw$tpr)
    auc = c(auc, roc.raw$auc)
    
    pvals = cbind(pvals, results$pvals)
    pvals.de = cbind(pvals.de, results$pvals[sim$de == 1])
    pvals.ee = cbind(pvals.ee, results$pvals[sim$de == 0])
  }
  
  update.sim(sim, tst.setting = 2, roc, auc, 
             pvals, pvals.de, pvals.ee)
}

QuasiSeq.tests = function(sim){
  print("Doing QuasiSeq tests.")
  
  roc.ql <- roc.qlshrink <- roc.qlspline <- NULL
  auc.ql <- auc.qlshrink <- auc.qlspline <- NULL
  pvals.ql <- pvals.qlshrink <- pvals.qlspline <- NULL
  pvals.de.ql <- pvals.de.qlshrink <- pvals.de.qlspline <- NULL
  pvals.ee.ql <- pvals.ee.qlshrink <- pvals.ee.qlspline <- NULL
  
  for(i in 1:dim(sim$phi)[2]){
    print(paste("  Testing with dispersions", colnames(sim$phi)[i]))
    design.list = vector("list", 2)
    design.list[[1]] = model.matrix(~sim$group)
    design.list[[2]] = rep(1, dim(sim$counts)[2])
    blahblahSPLINEblahblah = capture.output(
      fit <- suppressWarnings(QL.fit(counts = sim$counts, design.list = design.list, 
                                     Model = "NegBin", NBdisp = sim$phi[,i],
                                     log.offset = log(sim$norm),
                                     print.progress = F)))
    
    blahblahSPLINEblahblah = capture.output(
      res <- suppressWarnings(QL.results(fit, Plot = F)))
    
    allpv = res$P.values
    results.ql = data.frame(pvals = as.vector(allpv[["QL"]]), 
                            de = sim$de)
    results.qlshrink = data.frame(pvals = as.vector(allpv[["QLShrink"]]), 
                                  de = sim$de)  
    results.qlspline = data.frame(pvals = as.vector(allpv[["QLSpline"]]), 
                                  de = sim$de)

    roc.raw = results2roc(results.ql)
    roc.ql = cbind(roc.ql, roc.raw$tpr)
    auc.ql = c(auc.ql, roc.raw$auc)
    pvals.ql = cbind(pvals.ql, results.ql$pvals)
    pvals.de.ql = cbind(pvals.de.ql, results.ql$pvals[sim$de == 1])
    pvals.ee.ql = cbind(pvals.ee.ql, results.ql$pvals[sim$de == 0])
    
    roc.raw = results2roc(results.qlshrink)
    roc.qlshrink = cbind(roc.qlshrink, roc.raw$tpr)
    auc.qlshrink = c(auc.qlshrink, roc.raw$auc)
    pvals.qlshrink = cbind(pvals.qlshrink, results.qlshrink$pvals)
    pvals.de.qlshrink = cbind(pvals.de.qlshrink, 
                              results.qlshrink$pvals[sim$de == 1])
    pvals.ee.qlshrink = cbind(pvals.ee.qlshrink, 
                              results.qlshrink$pvals[sim$de == 0])
    
    roc.raw = results2roc(results.qlspline)
    roc.qlspline = cbind(roc.qlspline, roc.raw$tpr)
    auc.qlspline = c(auc.qlspline, roc.raw$auc)
    pvals.qlspline = cbind(pvals.qlspline, results.qlspline$pvals)
    pvals.de.qlspline = cbind(pvals.de.qlspline, 
                              results.qlspline$pvals[sim$de == 1])
    pvals.ee.qlspline = cbind(pvals.ee.qlspline, 
                              results.qlspline$pvals[sim$de == 0])
  }
  
  sim = update.sim(sim, tst.setting = 3, roc.ql, auc.ql, 
                   pvals.ql, pvals.de.ql, pvals.ee.ql)

  sim = update.sim(sim, tst.setting = 4, 
                   roc.qlshrink, auc.qlshrink, 
                   pvals.qlshrink, pvals.de.qlshrink, 
                   pvals.ee.qlshrink)

  sim = update.sim(sim, tst.setting = 5, roc.qlspline, 
                   auc.qlspline, pvals.qlspline, 
                   pvals.de.qlspline, pvals.ee.qlspline)
  sim
}
    
getDispersions = function(sim){
  print("Getting dispersions.")
  
  err(sim <- QL.dispersions(sim, dsp.setting = "QL"), 
      sim, "QL")

  err(sim <- DSS.dispersions(sim, dsp.setting = "DSS", 
      trend = F), sim, "DSS")

  err(sim <- wqCML.dispersions(sim, common = T, 
      dsp.setting = "wqCML: Common"), sim, "wqCML: Common")

  err(sim <- wqCML.dispersions(sim, common = F,
      dsp.setting = "wqCML: Tagwise"), sim, "wqCML: Tagwise")  

  err(sim <- APL.dispersions(sim, dsp.setting = "APL: Common", 
      method = "Common"), sim, "APL: Common")

  err(sim <- APL.dispersions(sim, dsp.setting = "APL: Trended", 
      method = "Trended"), sim, "APL: Trended")

  err(sim <- APL.dispersions(sim, dsp.setting = "APL: Tagwise", 
      method = "Tagwise"), sim, "APL: Tagwise")

  err(sim <- DESeq.dispersions(sim, dsp.setting = "DESeq: Trended", 
      sharingMode = "fit-only", fitType = "local"), sim, "DESeq: Trended")

  err(sim <- DESeq.dispersions(sim, dsp.setting = "DESeq: Maximum", 
      sharingMode = "maximum", fitType = "local"), sim, "DESeq: Maximum")
  
  err(sim <- DESeq.dispersions(sim, dsp.setting = "DESeq: None", 
      sharingMode = "gene-est-only"), sim, "DESeq: None")

  sim
}

doTests = function(sim){
  err(sim <- edgeR.test(sim), sim, "edgeR tests")
  err(sim <- DESeq.test(sim), sim, "DESeq tests")
  err(sim <- QuasiSeq.tests(sim), sim, "QuasiSeq tests")
  sim
}

oneSim = function(arg, nct){
  sim = sim.counts(arg)
  sim = getDispersions(sim)
  sim = doTests(sim)
  write(paste(as.roman(arg$sim.setting), nct, ":", 
              arg$nreps, "rep done"), 
        "progress.txt", append = T)
  print(paste(as.roman(arg$sim.setting), nct, ":", 
              arg$nreps, "rep done"))
  sim
}

oneSimSetting = function(arg){
  print(paste("Beginning new simulation setting:", 
              as.roman(arg$sim.setting)))
  write(paste(date()), "progress.txt", append = T)
  write(paste(as.roman(arg$sim.setting), "begun"), 
        "progress.txt", append = T)
  
  sims = mclapply(1:arg$nreps, 
                  function(nct){
                    err(oneSim(arg, nct), 
                               sim = arg,
                               info = paste("oneSim() rep", nct))},
                  mc.cores = arg$ncores)

  write(paste(as.roman(arg$sim.setting), "done\n"), 
        "progress.txt", append = T)
  print(paste("Done with simulation setting:", 
              as.roman(arg$sim.setting)))
  
  ret = list()
  new.names = paste("rep", 1:length(sims), sep="")
  
  counts = abind(lapply(sims, function(sim){sim$counts}), 
                 along = 3, new.names = new.names)
  saveRDS(counts, paste("../rds/counts", as.roman(arg$sim.setting),
                        ".rds", sep=""))
  
  delta = abind(lapply(sims, function(sim){sim$delta}), 
                along = 2, new.names = new.names)
  ret$delta = delta
  
  de = abind(lapply(sims, function(sim){sim$de}),
             along = 2, new.names = new.names)
  ret$de = de

  phi = abind(lapply(sims, function(sim){sim$phi}), 
              along = 3, new.names = new.names)
  phi = aperm(phi, c(2,1,3))
  ret$phi = phi
  
  auc = abind(lapply(sims, function(sim){sim$auc}), 
              along = 3, new.names = new.names)
  auc = aperm(auc, c(2, 1, 3))
  ret$auc = auc
  
  roc = abind(lapply(sims, function(sim){sim$roc}),
              along = 4, new.names = new.names)
  dimnames(roc)[[1]] = paste("tpr", 1:1000, sep = "")
  saveRDS(roc, paste("../rds/roc", as.roman(arg$sim.setting),
                        ".rds", sep=""))
  
  pvals = abind(lapply(sims, function(sim){sim$pvals}),
                along = 4, new.names = new.names)
  saveRDS(pvals, paste("../rds/pvals", as.roman(arg$sim.setting),
                        ".rds", sep=""))
  ret
}

runSimulations = function(){
  print("Running simulations.")
  args = readRDS("../rds/args.rds")
  simSettings = lapply(args, oneSimSetting)
  
  sim.names = as.character(as.roman(1:length(simSettings)))
  phi = abind(lapply(simSettings, function(x){x$phi}),
              along = 4, new.names = sim.names)
  phi = aperm(phi, c(4,1,2,3))
  saveRDS(phi, "../rds/phi.rds")
  
  auc = abind(lapply(simSettings, function(x){x$auc}),
              along = 4, new.names = sim.names)
  auc = aperm(auc, c(4,1,2,3))
  saveRDS(auc, "../rds/auc.rds")
  
  delta = abind(lapply(simSettings, function(x){x$delta}), 
                along = 3, new.names = sim.names)
  delta = aperm(delta, c(3,1,2))
  saveRDS(delta, "../rds/delta.rds")
  
  de = abind(lapply(simSettings, function(x){x$de}),
             along = 3, new.names = sim.names)
  de = aperm(de, c(3,1,2))
  saveRDS(de, "../rds/de.rds")
}

make.myframe = function(phi){
  pickrell_means = readRDS("../rds/true_meansI.rds")
  hammer_means = readRDS("../rds/true_meansIV.rds")
  pickrell_disps = readRDS("../rds/true_dispsI.rds")
  hammer_disps = readRDS("../rds/true_dispsIV.rds")
  
  pseudo_pk_means = geom.rowMeans(readRDS("../rds/countsI.rds")[,,1])
  pseudo_hm_means = geom.rowMeans(readRDS("../rds/countsIV.rds")[,,1])
  
  pseudo_pk_disps = phi[1,"True",,1]
  pseudo_hm_disps = phi[4,"True",,1]
  
  df = data.frame(
    ldisp = c(log(hammer_disps), log(pseudo_hm_disps), log(pickrell_disps), log(pseudo_pk_disps)),
    lmean = c(log(hammer_means), log(pseudo_hm_means), log(pickrell_means), log(pseudo_pk_means)))
  
  df$Dataset = c(rep("Hammer", length(hammer_disps) + length(pseudo_hm_disps)),
                 rep("Pickrell", length(pickrell_disps) + length(pseudo_pk_disps)))
  
  df$Status = c(rep("Actual Dataset", length(hammer_disps)),
                rep("Example Pseudo-dataset", length(pseudo_hm_disps)),
                rep("Actual Dataset", length(pickrell_disps)),
                rep("Example Pseudo-dataset", length(pseudo_pk_disps)))
  df
}

make.hists = function(df){  
  mdf = melt(df, id.vars = c("Dataset", "Status"))
  v = as.vector(mdf$variable)
  v[v == "ldisp"] = "Log dispersions"
  v[v == "lmean"] = "Log geometric mean counts"
  mdf$variable = v
  
  pl = ggplot(mdf, aes(x = value, y = ..density..)) + 
    geom_histogram(data = mdf, aes(fill=Dataset), binwidth=0.5, alpha=0.5, 
                   position="identity") + 
    facet_grid(Status ~ variable, scales = "free_x") + 
    xlab("") + ylab("") +
    scale_fill_manual(values = c("blue", "black")) + # SEPT 2013 REVISION 
    theme(text = element_text(family = "Helvetica", size = 12, colour= "black"),
          title= element_text(family = "Helvetica", size = 12),
          axis.text = element_text(family = "Helvetica", colour = 'black'),
          axis.line=element_line(),
          legend.title= element_text(family = "Helvetica", size=12, face="plain"))  
  
  
  ggsave("../fig/hists_slides.pdf", pl, width = 6, height = 6)
  
  pl = pl + theme(panel.background = element_rect(fill='white'),
                  panel.grid.major = element_line(color="lightgray"),
                  panel.border = element_rect(color="black", fill = NA))
  
  ggsave("../fig/hists.pdf", pl, width = 6, height = 6, dpi = 1600)
  
  pl = pl + theme(legend.position = "none")
  
  ggsave("../fig/hists_nolegend.pdf", pl, width = 6, height = 6, dpi = 1600)
}

meandisp.relation = function(df){
  pl = ggplot(df, aes(x = lmean, y = ldisp), na.rm=T) + 
    stat_binhex(na.rm = T, bins = 100, aes(fill=log(..count..))) +
    labs(fill="Log frequency") +
    facet_grid(Status~Dataset) + 
    xlab("Log geometric mean counts") + 
    ylab("Log dispersions") + 
    scale_fill_continuous(low = "grey80", high = "black") +
    theme(text = element_text(family = "Helvetica", size = 12, colour= "black"),
          axis.text = element_text(family = "Helvetica", colour = 'black'),
          axis.line=element_line(),
          legend.title= element_text(family = "Helvetica", size=12, face="plain"))
  
  ggsave("../fig/meandisp_scatter_slides.pdf", pl, 
         width = 8, height = 8, dpi = 1600)
  
  pl = pl + theme(panel.background = element_rect(fill='white'),
                  panel.grid.major = element_line(color="lightgray"),
                  panel.border = element_rect(color="black", fill = NA))
  
  ggsave("../fig/meandisp_scatter.pdf", pl, 
         width = 8, height = 8, dpi = 1600)
  
  pl = pl + theme(legend.position = "none")
  
  ggsave("../fig/meandisp_scatter_nolegend.pdf", pl, 
         width = 8, height = 8, dpi = 1600)
}

example.roc = function(){
  rocI = readRDS("../rds/rocI.rds")[,"QL", "test1","rep1"]
  rocVI = readRDS("../rds/rocVI.rds")[,"True", "test5", "rep1"]
  nms = c("Simulation setting I, QL dispersions, edgeR test", 
          "Simulation setting VI, True dispersions, QLSpline test")
  df = data.frame(TPR = c(0, rocI, 0, rocVI),
                  FPR = rep(0:1000/1000, 2),
                  Example = c(rep(nms, each = 1001)))

  pl = ggplot(df, aes(x = FPR, y = TPR), na.rm=TRUE) + 
    geom_line(aes(color = factor(Example))) +
    xlab("False positive rate (FPR)") + 
    ylab("True positive rate (TPR)") +
    labs(color = "Example ROC curves") + 
    geom_vline(xintercept = 0.1, linetype = "dotted")+
    theme(text = element_text(family = "Helvetica", size = 12, colour= "black", face="plain"),
          legend.title= element_text(family = "Helvetica", size=12, face="plain"),
          axis.text = element_text(family = "Helvetica", colour = 'black'),
          axis.line=element_line(),
          legend.position="top",
          legend.direction="vertical")
  
  ggsave("../fig/roc_slides.pdf", pl, width = 5, height = 5, dpi = 1600)
  
  pl = pl + theme(panel.background = element_rect(fill='white'),
                  panel.grid.major = element_line(color="lightgray"),
                  panel.border = element_rect(color="black", fill = NA))
  
  ggsave("../fig/roc.pdf", pl, width = 5, height = 5, dpi = 1600)
  
  pl = pl + theme(legend.position = "none")
  
  ggsave("../fig/roc_nolegend.pdf", pl, width = 5, height = 5, dpi = 1600)
}

plot.err = function(phmsem, filename, filename.slides, ylab){
  dimnames(phmsem)[[2]][dimnames(phmsem)[[2]] == "QL"] = "QL: None"
  dimnames(phmsem)[[2]][dimnames(phmsem)[[2]] == "DSS"] = "DSS: Tagwise"
  
  phmsef = adply(phmsem, 1:length(dim(phmsem)))
  colnames(phmsef) = c("sim", "dsp", "rep", "mse")
  phmsef$package = factor(rep(gsub(": .*", "", dimnames(phmsem)[[2]]) , 
                          each = length(dimnames(phmsem)[[1]])),
                          levels = c("True", "QL", "DSS", "wqCML", "APL", "DESeq"))
  
  pl = ggplot(phmsef, aes(x = dsp, y = mse), na.rm=TRUE) + 
         geom_boxplot(aes(x = dsp, y = mse)) + 
         facet_grid(sim ~ package, scales = "free_x", space = "free_x") + 
         xlab("\nDispersion shrinkage type") + ylab(ylab) +
         theme(axis.text.x = element_text(family = "Helvetica", angle = -45, hjust = 0)) + 
         scale_x_discrete(labels = function(x){gsub(".*: ", "", x)})+
    theme(text = element_text(family = "Helvetica", size = 12, colour= "black"),
          axis.text = element_text(family = "Helvetica", colour = 'black'),
          axis.line=element_line())
  
  ggsave(paste("../fig/",filename.slides, sep=""), pl, width = 8, height = 8, dpi = 1600)
  
  pl = pl + theme(panel.background = element_rect(fill='white'),
                  panel.grid.major = element_line(color="lightgray"),
                  panel.border = element_rect(color="black", fill = NA))
  
  ggsave(paste("../fig/", filename, sep=""), pl, width = 8, height = 8, dpi = 1600)
}

phi.vs.phi = function(phi, sim.setting = 1){
  df = melt(log(phi[sim.setting,,,1]))
  colnames(df) = c("Method", "Gene", "Estimated")
  
  true = df[df$Method == "True",]$Estimated
  
  df = df[df$Method != "True",]
  df = df[df$Method != "APL: Common",]
  df = df[df$Method != "wqCML: Common",]
  
  df$True = rep(true, each = 8)
  
  lmean = log(geom.rowMeans(readRDS(paste("../rds/counts",
                                          as.character(as.roman(sim.setting)),
                                          ".rds", sep=""))[,,1]))

  lmean = cut(lmean, breaks=quantile(lmean, probs = c(0, .5, 1)), 
              include.lowest = T) 

  df$lmean = rep(lmean, each = 8)
  
  df = df[df$True != min(df$True),]
  
  lim_x = quantile(df$True, probs = c(0, 1))
  lim_y = quantile(df$Estimated, probs = c(0,1))
  lim = c(-10, 1)
  
  pl = ggplot(df, aes(x = True, y = Estimated), na.rm=TRUE) + 
        stat_binhex(na.rm=TRUE, bins=100, 

#                   aes(fill=lmean, alpha = log(..count..))) +  # JULY 2013 SUBMISSION
#       scale_alpha(range = c(0.25, 1)) + # JULY 2013 SUBMISSION
        
                    aes(fill=lmean), alpha = 0.5) + # SEPT 2013 REVISION
        scale_fill_manual(values = c("black", "blue")) + # SEPT 2013 REVISION 

        xlim(lim_x) + ylim(lim_x) + 
        facet_wrap(~Method) + xlab("\nTrue log dispersions") + 
        ylab("Estimated log dispersions\n") + 
        labs(fill = "Log geometric mean count") + 
        geom_abline() +  
        guides(fill = guide_legend(order = 1)

               ) + # SEPT 2013 REVISION
            #   , alpha = guide_legend(order = 2)) + # JULY 2013 SUBMISSION

    theme(text = element_text(family = "Helvetica", size = 12, colour= "black"),
          axis.text = element_text(family = "Helvetica", colour = 'black'),
          axis.line=element_line(),
          legend.title= element_text(family = "Helvetica", size=12, face="plain"))
  
  ggsave(paste("../fig/phi_vs_phi_", 
               as.character(as.roman(sim.setting)),
               "_slides.pdf", sep=""), pl, width = 8, height = 8, dpi = 1600)
  
  pl = pl + theme(panel.background = element_rect(fill='white'),
                  panel.grid.major = element_line(color="lightgray"),
                  panel.border = element_rect(color="black", fill = NA))
  
  ggsave(paste("../fig/phi_vs_phi_", 
                 as.character(as.roman(sim.setting)),
                 ".pdf", sep=""), pl, width = 8, height = 8, dpi = 1600)
  
  pl = pl + theme(legend.position = "none")
  
  ggsave(paste("../fig/phi_vs_phi_", 
               as.character(as.roman(sim.setting)),
               "_nolegend.pdf", sep=""), pl, width = 8, height = 8, dpi = 1600)
}

plot.auc = function(){
  aucm = readRDS("../rds/auc.rds")
  dimnames(aucm)[[2]][dimnames(aucm)[[2]] == "True"] = "True: None"
  dimnames(aucm)[[2]][dimnames(aucm)[[2]] == "QL"] = "QL: None"
  dimnames(aucm)[[2]][dimnames(aucm)[[2]] == "DSS"] = "DSS: Tagwise"
  dimnames(aucm)[[3]] = c("edgeR", "DESeq", "QL", "QLShrink", "QLSpline")
  
  aucf = adply(aucm, 1:length(dim(aucm)))
  colnames(aucf) = c("sim", "dsp", "test", "rep", "auc")
  aucf$package = factor(rep(gsub(": .*", "", dimnames(aucm)[[2]]) , 
                        each = length(dimnames(aucm)[[1]])),
                        levels = c("True", "QL", "DSS", "wqCML", "APL", "DESeq"))
  
  filepath = "../fig/auc"
  for(i in 1:dim(aucm)[1]){
    df = aucf[aucf$sim == as.character(as.roman(i)),]
    filename = paste(filepath, i, ".pdf", sep="")
    
    pl = ggplot(df, aes(x = dsp, y = auc), na.rm=TRUE) + 
           geom_boxplot(aes(x = dsp, y = auc)) + 
           facet_grid(test ~ package, scales = "free_x", space = "free_x") + 
           xlab("\nDispersion shrinkage type") + ylab("Area under Receiver Operating Characteristic (ROC) curve \n(false positive rate < 0.1)\n") +
           scale_x_discrete(labels = function(x){gsub(".*: ", "", x)})+
      theme(axis.text.x = element_text(family = "Helvetica", angle = -45, hjust = 0),
            text = element_text(family = "Helvetica", size = 12, colour= "black"),
            axis.text = element_text(family = "Helvetica", colour = 'black'),
            axis.line=element_line(),
            legend.position="none")
      
    ggsave(paste(filepath, i, "_slides.pdf", sep=""),
           pl, width = 8, height = 8, dpi = 1600)
    
    pl = pl + theme(panel.background = element_rect(fill='white'),
                    panel.grid.major = element_line(color="lightgray"),
                    panel.border = element_rect(color="black", fill = NA)) +
              ylab("Area under ROC curve (false positive rate < 0.1)\n")
    
    ggsave(filename, pl, width = 8, height = 8, dpi = 1600)
  }
}

point.est = function(){
  phi = readRDS("../rds/phi.rds")
  df = make.myframe(phi)
  make.hists(df)
  meandisp.relation(df)
  
  se = sweep(phi[,-1,,]/(phi[,-1,,] + 1), 
                 MARGIN = c(1,3,4), 
                 STATS = phi[,1,,]/(phi[,1,,] + 1),
                 FUN = "-")^2
  
  mse = apply(se, c(1,2,4), mean)
    
  plot.err(mse, "mse.pdf", "mse_slides.pdf", "Mean squared error\n")
  
  for(i in 1:dim(phi)[1])
    phi.vs.phi(phi, i)
}

print.results = function(){
  print("Printing results.")
  example.roc()
  point.est()
  plot.auc()
  system("./plosone_figures.sh")
}

workflow = function(){
  initialize()
  computeInput()
  runSimulations()
  print.results()
  warnings()
}