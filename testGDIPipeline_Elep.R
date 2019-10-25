# test gdi pileine with lepidum data

library(gdiPipeline)
library(ape)
library(geiger)
library(stringr)
library(phytools)

# pipeline parameters

wd="/Users/dmacguigan/Documents/NearLab/LepidumProject/BPP/GDI"
treefile="Elep.tree"
map="Elep_allLoci.Imap.txt"
priors="priors.txt"
heredity = "heredity.txt"
loci = "Elep_allLoci.txt"
ctl = "ctlTemplate.ctl"
plotColors = c(brewer.pal(12, "Paired"))

# pipeline steps

BPPCtlTemplate(wd)

bppInputs(wd, treefile, map,
          priors,
          heredity, loci,
          ctl, 14)

bppTaskFile(wd)

gdi <- bppSummarizeGDI(wd, plotColors)
