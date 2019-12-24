# gdiPipeline
## IN DEVELOPMENT

R package to run GDI species delimitation analyses using BPP

To download and install this package
```
install.packages("devtools")
library(devtools)
install_github("dmacguigan/gdiPipeline")
library(gdiPipeline)
```

Quick tutorial
```
# test gdi pileine with E. lepidum data

library(gdiPipeline)

# pipeline parameters

wd="/Users/dmacguigan/Documents/NearLab/LepidumProject/BPP/GDI/bpp-4.1_reps_2"
treefile="Elep.tree"
map="Elep_allLoci.Imap.txt"
priors="priors.txt"
heredity = "heredity.txt"
loci = "Elep_allLoci.txt"
ctl = "ctlTemplate.ctl"
plotColors = c(brewer.pal(12, "Paired"), "black")
nLoci = 14
threads = 10
nreps = 4

# pipeline steps

# create BPP control file template for bppInputs function
BPPCtlTemplate(wd)

# create directories and files for BPP runs
# specify priors (txt file, each line with a specific combination)
# specify maximally split tree
bppInputs(wd, treefile, map,
          priors,
          heredity, loci,
          ctl, nLoci,
          threads, nreps)

# create task file, each line with a different command to run BPP
# use this task file in conjunction with job array on computing cluster (e.g. https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/dsq/)
bppTaskFile(wd)

# summarize the bpp results
# caluclate GDI
# plot GDI for each replicate run and binning replicate runs
gdi <- bppSummarizeGDI(wd, plotColors, nreps)

# plot GDI estimates for all species in one figure
plotByPrior(gdi, wd, nreps, priors, plotWidth=16, plotHeight=16)
```
