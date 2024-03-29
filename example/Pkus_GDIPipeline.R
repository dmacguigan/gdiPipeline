# GDI pipeline with Percina kusha example data

library(gdiPipeline)
library(RColorBrewer)

############################################################################################################################################
# pipeline parameters

wd="G:/NearLab/Pkus/analyses/BPP/m100p_BPPSetup"
treefile="Pkus.tree" # make sure newick tree has semicolon at the line end
map="Pkus.Imap.txt"
priors="priors.txt"
heredity = "heredity.txt"
loci = "m100p_10loci.phy.txt"
ctl = "ctlTemplate.ctl"
plotColors = c(brewer.pal(12, "Paired"), "black")
nLoci = 10
threads = 8
nreps = 4

############################################################################################################################################
# pipeline steps

############################################################################################################################################
# STEP 1
# create BPP control file template for bppInputs function
BPPCtlTemplate(wd)

############################################################################################################################################
# STEP 2
# modify burnin, samplefreq, or nsample parameters in the ctlTemplate.ctl file if necessary before proceeding
# may also want to specify finetune values based on preliminary runs if mixing/convergence is an issue

############################################################################################################################################
# STEP 3
# create directories and files for BPP runs
# specify priors (txt file, each line with a specific combination)
# specify maximally split tree
bppInputs(wd, treefile, map,
          priors,
          heredity, loci,
          ctl, nLoci,
          threads, nreps)

############################################################################################################################################
# STEP 4
# create task file, each line with a different command to run BPP
# use this task file in conjunction with job array on computing cluster (e.g. https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)
bppTaskFile(wd)

############################################################################################################################################
# STEP 5
# run the task file generated by bppTaskFile
# ideally parallelize using a job array setup (e.g. https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)

############################################################################################################################################
# STEP 6
# assess mixing and convergence of BPP analyses
# diagnostic files are written in each model directory
# this may take a while to run depending on the number of analyses and the MCMC chain lengths
checkConvergence(wd=wd, nreps=nreps)

############################################################################################################################################
# STEP 7
# summarize the bpp results
# caluclate GDI
# plot GDI for each replicate run and binning replicate runs
gdi <- bppSummarizeGDI(wd, plotColors, nreps)

# plot GDI estimates for all species in one figure
plotByPrior(gdi, wd, nreps, priors, plotWidth=8, plotHeight=5)
