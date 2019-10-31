# test gdi pileine with lepidum data

library(gdiPipeline)
library(RColorBrewer)
library(geiger)
library(stringr)
library(phytools)

# pipeline parameters

wd="/Users/dmacguigan/Documents/NearLab/LepidumProject/BPP/GDI/bpp-4.1_reps"
treefile="Elep.tree"
map="Elep_allLoci.Imap.txt"
priors="priors.txt"
heredity = "heredity.txt"
loci = "Elep_allLoci.txt"
ctl = "ctlTemplate.ctl"
plotColors = c(brewer.pal(12, "Paired"))
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

gdi <- bppSummarizeGDI(wd, plotColors)
