# create inputs for BPP analyses

BPPCtlTemplate <- function(wd){
  setwd(wd)
  fileConn<-file("ctlTemplate.ctl", "wb")
  writeLines(c("  seed =  -1",
  "  seqfile = loci.txt",
  "  Imapfile = Imap.txt",
  "  outfile = out.txt",
  "  mcmcfile = mcmc.txt",
  "  speciesdelimitation = 0",
  "  speciestree = 0",
  "  speciesmodelprior = 1",
  "  species&tree = ",
  "                 nsamples",
  "                 phy;",
  "  diploid = 1",
  "  cleandata = 0",
  "  usedata = 1",
  "  nloci = n",
  "  thetaprior = a b",
  "  tauprior = a b",
  "  heredity = 2 heredity.txt",
  "  locusrate = 1 2.0",
  "  finetune =  1: 5 0.001 0.001  0.001 0.3 0.33 1.0",
  "  threads = nThreads",
  "  print = 1 0 0 0",
  "  burnin = 10000",
  "  sampfreq = 50",
  "  nsample = 10000"), fileConn)
  close(fileConn)
}

bppInputs <- function(wd, treefile, map, priors, heredity, loci, ctl, nloci, threads) {
  # set working directory
  setwd(wd)

  # read in priors
  prior_df <- read.table(priors, header=FALSE)

  # read in the tree
  tree <- read.tree(treefile)

  #get tip labels
  tipLabs <- (tree$tip.label)

  # get internal node numbers, skipping root node
  n <- length(tipLabs) + tree$Nnode
  n <- (length(tipLabs) + 2):n

  # first create directories for fully split model (always model 1)
  # for each combination of priors
  dir.create("model1")
  for(j in 1:nrow(prior_df)){
    newDir = paste("model1/tau",paste(prior_df[j,1:2], collapse="-"),"theta",paste(prior_df[j,3:4], collapse="-"), sep="")
    dir.create(newDir)
    file.copy(loci, paste(newDir,"/loci.txt", sep=""))
    file.copy(heredity, paste(newDir,"/heredity.txt", sep=""))
    file.copy(map, paste(newDir,"/Imap.txt", sep=""))
    file.copy(ctl, paste(newDir,"/bpp.ctl", sep=""))
    # write the tree for this delimitation model to a file
    write.tree(tree, file="model1.tree")
    tree_newick <- readLines("model1.tree")
    # modify the BPP control file and the map file
    ctlTxt <- readLines(paste(newDir,"/bpp.ctl",sep=""))
    map_df <- read.table(paste(newDir,"/Imap.txt", sep=""))
    map_table <- table(map_df[,2])
    taxa_names <- paste(names(map_table), collapse=" ")
    taxa_counts <- paste(map_table, collapse=" ")
    taxa <- (tree$tip.label)
    ctlTxt <- sub("phy;", tree_newick, ctlTxt)
    ctlTxt <- sub("species&tree = ", paste("species&tree = ", length(tipLabs), " ", taxa_names, sep=""), ctlTxt)
    ctlTxt <- sub("nsamples", taxa_counts, ctlTxt)
    ctlTxt <- sub("nThreads", threads, ctlTxt)
    ctlTxt <- sub("diploid = 1", paste("diploid = ", paste(rep(1, length(tipLabs)), collapse=" "), sep = ""), ctlTxt)
    ctlTxt <- sub("thetaprior = a b", paste("thetaprior = ", paste(prior_df[j,3:4], collapse=" "), " E", sep=""), ctlTxt)
    ctlTxt <- sub("tauprior = a b", paste("tauprior = ", paste(prior_df[j,1:2], collapse=" "), sep=""), ctlTxt)
    ctlTxt <- sub("nloci = n", paste("nloci = ", nloci, sep=""), ctlTxt)
    f <- file(paste(newDir,"/bpp.ctl",sep=""), "wb")
    writeLines(ctlTxt, f)
    tbl <- file(paste(newDir,"/Imap.txt", sep=""), "wb")
    write.table(map_df, file=tbl, row.names = FALSE, col.names = FALSE, quote=FALSE)
  }

  # for each clade in the tree
  # first
  count = 2
  for(i in n){
    # create new directory for each sp delim mode
    modelDir <- paste("model", count, sep="")
    dir.create(modelDir)
    count = count + 1
    # get tips for clade i
    t <- tips(tree, node=i)
    t_prune <- t[2:length(t)]
    t_keep <- t[1]
    # create new name for collpased node
    newName <- paste(t, collapse="")
    # drop all but one tip in the clade
    tempT <- drop.tip(tree, t_prune)
    # rename the remaining tip
    tempT$tip.label[tempT$tip.label==t_keep] <- newName
    # write new tree to a file
    write.tree(tempT, file=paste(modelDir, ".tree", sep=""))
    # for each combination of priors
    for(j in 1:nrow(prior_df)){
      newDir = paste(modelDir, "/tau",paste(prior_df[j,1:2], collapse="-"),"theta",paste(prior_df[j,3:4], collapse="-"), sep="")
      dir.create(newDir)
      file.copy(loci, paste(newDir,"/loci.txt", sep=""))
      file.copy(heredity, paste(newDir,"/heredity.txt", sep=""))
      file.copy(map, paste(newDir,"/Imap.txt", sep=""))
      file.copy(ctl, paste(newDir,"/bpp.ctl", sep=""))
      # write the tree for this delimitation model to a file
      tree_newick <- readLines(paste(modelDir, ".tree", sep=""))
      # read in new tree
      newTree <- read.tree(paste(modelDir, ".tree", sep=""))
      newTreeTaxa <- (newTree$tip.label)
      # modify the BPP control file and the map file
      ctlTxt <- readLines(paste(newDir,"/bpp.ctl",sep=""))
      map_df <- read.table(paste(newDir,"/Imap.txt", sep=""))
      # replace names in the map file
      for(k in t){
        map_df[2] <- lapply(map_df[2], gsub, pattern = k, replacement = "fakeName")
      }
      map_df[2] <- lapply(map_df[2], gsub, pattern = "fakeName", replacement = newName)
      map_table <- table(map_df[,2])
      taxa_names <- paste(names(map_table), collapse=" ")
      taxa_counts <- paste(map_table, collapse=" ")
      taxa <- (tree$tip.label)
      ctlTxt <- sub("phy;", tree_newick, ctlTxt)
      ctlTxt <- sub("species&tree = ", paste("species&tree = ", length(map_table), " ", taxa_names, sep=""), ctlTxt)
      ctlTxt <- sub("nsamples", taxa_counts, ctlTxt)
      ctlTxt <- sub("thetaprior = a b", paste("thetaprior = ", paste(prior_df[j,3:4], collapse=" "), " e", sep=""), ctlTxt)
      ctlTxt <- sub("diploid = 1", paste("diploid = ", paste(rep(1, length(newTreeTaxa)), collapse=" "), sep = ""), ctlTxt)
      ctlTxt <- sub("tauprior = a b", paste("tauprior = ", paste(prior_df[j,1:2], collapse=" "), sep=""), ctlTxt)
      ctlTxt <- sub("nloci = n", paste("nloci = ", nloci, sep=""), ctlTxt)
      ctlTxt <- sub("nThreads", threads, ctlTxt)
      writeLines(ctlTxt, paste(newDir,"/bpp.ctl",sep=""))
      write.table(map_df, file=paste(newDir,"/Imap.txt", sep=""), row.names = FALSE, col.names = FALSE, quote=FALSE)
    }
  }
}

bppTaskFile <- function(wd, col) {
  dirs <- list.dirs(path=wd)
  dirs <- grep("tau", dirs, value=TRUE)
  dirs <- gsub("/.*/model", "./model", dirs)
  commands <- NULL
  for(i in dirs){
    newCommand = paste("cd ", i, ";bpp --cfile bpp.ctl", sep="")
    commands = c(commands, newCommand)
  }
  fileConn<-file("BPPTaskFile.txt")
  writeLines(commands, con = fileConn)
  close(fileConn)
}

bppSummarizeGDI <- function(wd, col) {
  setwd(wd)
  # get all BPP dirs in the wd
  dirs <- list.dirs(path=wd)
  dirs <- grep("tau", dirs, value=TRUE)
  # list to store all gdi estimates
  allGDIList <- vector(mode = "list", length = length(dirs))
  count=1
  for(i in dirs) {
    # read in the posterior
    mcmc <- read.table(paste(i, "/mcmc.txt", sep=""), header=TRUE)
    p <- strsplit(i, "/")
    model <- p[[1]][length(p[[1]])-1]
    # read in the model tree
    tree <- read.tree(paste(model, ".tree", sep=""))
    taxa <- tree$tip.label
    thetaList <- vector(mode = "list", length = length(taxa))
    tauList <- vector(mode = "list", length = length(taxa))
    gdiList <- vector(mode = "list", length = length(taxa))
    gdiTaxa <- vector(mode = "list", length = length(taxa))
    cols <- colnames(mcmc)
    # for each taxon in the tree
    for(j in 1:length(taxa)){
      # get theta that matches taxon
      thetaList[[j]] <- mcmc[,j+1]
      temp <- sapply(taxa, function(x) (grepl(x, names(mcmc)[j+1], fixed = TRUE)))
      sp <- taxa[temp]
      # get sister species or node of taxon
      sis <- getSisters(tree, sp, mode="label")
      if(is.character(sis[[1]])){
        node <- findMRCA(tree, tips=c(sis[[1]],sp), type="node")
        index <- ncol(mcmc) - ((length(taxa) + tree$Nnode) + 1 - node)
        tauList[[j]] <- mcmc[,index]
      } else {
        node <- sis[[1]] - 1
        index <- ncol(mcmc) - ((length(taxa) + tree$Nnode) + 1 - node)
        tauList[[j]] <- mcmc[,index]
      }
      # calculate gdi
      gdiList[[j]] <- (1-exp((-2*tauList[[j]]/thetaList[[j]])))
      # save list of species
      gdiTaxa[[j]] <- sp
    }
    # plot gdi
    p <- strsplit(i, "/")
    priors <- p[[1]][length(p[[1]])]
    print(model)
    print(priors)
    dat <- data.frame(matrix(unlist(gdiList), ncol=length(gdiList), byrow=F))
    colnames(dat) <- unlist(gdiTaxa)
    pdf(file=paste(model, "_", priors, "_gdi.pdf", sep=""), width=8, height = 6)
    p <- ggplot(data = melt(dat), aes(x=value, color=variable, fill=variable)) +
      xlim(c(0,1)) +
      geom_vline(xintercept = 0.2, lty=2)+
      geom_vline(xintercept = 0.7, lty=2)+
      geom_density(alpha=0.1) +
      scale_color_manual(values=col) +
      scale_fill_manual(values=col) +
      labs(color="species", fill="species", x="GDI")+
      ggtitle(label=paste(model, priors, sep=" ")) +
      annotate("text", label="species", x=0.85, y=Inf, hjust=0.5, vjust=1) +
      annotate("text", label="ambiguous", x=0.45, y=Inf, hjust=0.5, vjust=1) +
      annotate("text", label="populations", x=0.1, y=Inf, hjust=0.5, vjust=1) +
      theme_minimal()
    print(p)
    dev.off()
    allGDIList[[count]] <- dat
    names(allGDIList)[count] <- paste(model, priors, sep="")
    count=count+1
  }
  return(allGDIList)
}
