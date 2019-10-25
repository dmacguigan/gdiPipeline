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
# test gdi pileine with lepidum data

library(gdiPipeline)

# pipeline parameters

wd="/Users/dmacguigan/Documents/NearLab/LepidumProject/BPP/GDI"
treefile="Elep.tree"
map="Elep_allLoci.Imap.txt"
priors="priors.txt"
heredity = "heredity.txt"
loci = "Elep_allLoci.txt"
ctl = "ctlTemplate.ctl" # this file will be created below by the function BPPCtlTemplate
plotColors = c(brewer.pal(12, "Paired"))
nloci = 14

# pipeline steps

BPPCtlTemplate(wd)

bppInputs(wd, treefile, map,
          priors,
          heredity, loci,
          ctl, nloci)

bppTaskFile(wd)

gdi <- bppSummarizeGDI(wd, plotColors)
```
