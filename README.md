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

BPPCtlTemplate(wd)

bppInputs(wd, treefile, map,
          priors,
          heredity, loci,
          ctl, nLoci,
          threads, nreps)

bppTaskFile(wd)

gdi <- bppSummarizeGDI(wd, plotColors, nreps)

plotByPrior(gdi, wd, nreps, priors, plotWidth=16, plotHeight=16)
```
