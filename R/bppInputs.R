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

bppInputs <- function(wd, treefile, map, priors, heredity, loci, ctl, nloci, threads, nreps) {
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
    close(f)
    close(tbl)
    for(z in 1:nreps){
      dir.create(paste(newDir, z, sep="-"))
      fileList <- list.files(newDir, full.names=TRUE)
      file.copy(fileList, paste(newDir, z, sep="-"))
    }
    # delete the first directory
    unlink(newDir, recursive = TRUE)
  }

  # for each clade in the tree
  # first
  count = 2
  # if there are more than 2 tips in the guide tree...
  if(length(tipLabs) > 2){
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
		  for(z in 1:nreps){
			dir.create(paste(newDir, z, sep="-"))
			fileList <- list.files(newDir, full.names=TRUE)
			file.copy(fileList, paste(newDir, z, sep="-"))
		  }
		  # delete the first directory
		  unlink(newDir, recursive = TRUE)
		}
	}
  }
}

bppTaskFile <- function(wd) {
  dirs <- list.dirs(path=wd, full.names=FALSE)
  dirs <- grep("tau", dirs, value=TRUE)
  commands <- NULL
  for(i in dirs){
    newCommand = paste("cd ", i, ";bpp --cfile bpp.ctl > ", gsub("model[0-9]+/", "", i), ".out 2> ", gsub("model[0-9]+/", "", i), ".error", sep="")
    commands = c(commands, newCommand)
  }
  fileConn<-file("BPPTaskFile.txt")
  writeLines(commands, con = fileConn)
  close(fileConn)
}

bppSummarizeGDI <- function(wd, col, nreps, burnin) {

  replicate = 0
  setwd(wd)

  # get all BPP dirs in the wd
  dirs <- list.dirs(path=wd)
  dirs <- grep("tau", dirs, value=TRUE)

  # list to store all gdi estimates
  allGDIList <- vector(mode = "list", length = length(dirs))
  count=1
  for(i in dirs) {
    replicate = replicate + 1
    rep <- strsplit(i, "/")
    rep <- rep[[1]][length(rep[[1]])]
    # read in population map from std out file
    con <- file(paste(i, "/", rep, ".out", sep=""), open="r")
    lines <- c()
    read=FALSE
    while(length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      if (startsWith(line, "Map of")) {
        read=TRUE
        next
      } else if (startsWith(line, "Generating")){
        break
      } else if (read == TRUE){
        lines <- c(lines, line)
      } else if ( length(line) == 0 ) {
        break
      }
    }
    close(con)

    # split lines and extract taxa names
    lines <- head(lines,-1)
    lines <- lines[2:length(lines)]
    lines <- strsplit(lines, " ")
    lines <- lapply(lines, function(x){x[!x ==""]})
    species <- unlist(lapply(lines, '[[', 2))

    # read in the posterior
    mcmc <- read.table(paste(i, "/mcmc.txt", sep=""), header=TRUE)
    prior <- strsplit(i, "/")
    model <- prior[[1]][length(prior[[1]])-1]
	
	# remove burnin as percentage of the chain
	gen <- nrow(mcmc)
	if(burnin > 0){
		burnin_gen <- round(gen*burnin, digits=0)
		mcmc <- mcmc[-c(1:burnin_gen),]
	}

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
      # get species index
      spIndex <- match(taxa[j], species)

      # get theta that matches taxon
      thetaList[[j]] <- mcmc[,spIndex+1]
      sp <- species[spIndex]

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
    priors <- strsplit(i, "/")
    priors <- priors[[1]][length(priors[[1]])]
    print(model)
    print(priors)
    dat_theta <- data.frame(matrix(unlist(thetaList), ncol=length(thetaList), byrow=F))
    dat_tau <- data.frame(matrix(unlist(tauList), ncol=length(tauList), byrow=F))
    dat <- data.frame(matrix(unlist(gdiList), ncol=length(gdiList), byrow=F))
    colnames(dat_theta) <- unlist(gdiTaxa)
    colnames(dat_tau) <- unlist(gdiTaxa)
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

    pdf(file=paste(model, "_", priors, "_gdi_boxplot.pdf", sep=""), width=8, height = 6)
    p <- ggplot(data = melt(dat), aes(y=value, x=variable, fill=variable)) +
      ylim(c(0,1)) +
      geom_hline(yintercept = 0.2, lty=2)+
      geom_hline(yintercept = 0.7, lty=2)+
      geom_boxplot() +
      scale_fill_manual(values=col) +
      labs(fill="species", y="GDI", x="species")+
      ggtitle(label=paste(model, priors, sep=" ")) +
      annotate("text", label="species", y=0.85, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      annotate("text", label="ambiguous", y=0.45, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      annotate("text", label="populations", y=0.1, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")
    print(p)
    dev.off()

    # append GDI data to list
    allGDIList[[count]] <- dat
    names(allGDIList)[count] <- paste(model, priors, sep="")

    if (replicate == nreps){
      replicate = 0
      for (i in 1:(nreps-1)){
        dat <- rbind(dat, allGDIList[[count-i]])
      }
      colnames(dat) <- unlist(gdiTaxa)
      priors <- substr(priors, 1, nchar(priors)-2)

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

      pdf(file=paste(model, "_", priors, "_gdi_boxplot.pdf", sep=""), width=8, height = 6)
      p <- ggplot(data = melt(dat), aes(y=value, x=variable, fill=variable)) +
        ylim(c(0,1)) +
        geom_hline(yintercept = 0.2, lty=2)+
        geom_hline(yintercept = 0.7, lty=2)+
        geom_boxplot() +
        scale_fill_manual(values=col) +
        labs(fill="species", y="GDI", x="species")+
        ggtitle(label=paste(model, priors, sep=" ")) +
        annotate("text", label="species", y=0.85, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
        annotate("text", label="ambiguous", y=0.45, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
        annotate("text", label="populations", y=0.1, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")
      print(p)
      dev.off()

    }
    count=count+1
  }
  return(allGDIList)

}

mergeReplicates <- function(wd, col, nreps, priors) {

}

plotByPrior <- function(gdiDat, wd, nreps, priors, plotWidth, plotHeight) {

  setwd(wd)	

  # read in priors
  prior_df <- read.table(priors, header=FALSE)
  n_priors <- nrow(prior_df)
  allGDIList_byPrior <- vector(mode = "list", length = n_priors)

  count = 0
  replicate = 0
  p_count = 0

  for(i in 1:length(gdiDat)) {
    replicate = replicate + 1
    count = count + 1
    if(replicate == nreps){
      p_count = p_count + 1
      replicate = 0
      dat = gdiDat[[i]]
      for (j in 1:(nreps-1)){
        dat <- rbind(dat, gdiDat[[count-j]])
      }
      # create data frame of GDI values for each prior, only include one instance of each taxon
      for (k in 1:ncol(dat)){
        taxa <- colnames(allGDIList_byPrior[[p_count]])
        if (!((colnames(dat))[k] %in% taxa)){
          allGDIList_byPrior[[p_count]] <- cbind(allGDIList_byPrior[[p_count]], dat[,k])
          colnames(allGDIList_byPrior[[p_count]])[ncol(allGDIList_byPrior[[p_count]])] <- colnames(dat)[k]
        }
      }
    }
    if(p_count == n_priors){
      p_count = 0
    }
  }

  for (i in 1:length(allGDIList_byPrior)) {
    dat = as.data.frame(allGDIList_byPrior[[i]])
    pdf(file=paste("priors-tau", prior_df[i,1], prior_df[i,2], "theta", prior_df[i,3], prior_df[i,4], "_boxplot.pdf", sep=""), width=plotWidth, height=plotHeight)
    p <- ggplot(data = melt(dat), aes(y=value, x=variable)) +
      ylim(c(0,1)) +
      geom_hline(yintercept = 0.2, lty=2)+
      geom_hline(yintercept = 0.7, lty=2)+
      geom_boxplot() +
      labs(y="GDI", x="species")+
      ggtitle(label=paste("priors: tau(", prior_df[i,1], ",", prior_df[i,2], "), theta(", prior_df[i,3], ",", prior_df[i,4], ")", sep="")) +
      annotate("text", label="species", y=0.85, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      annotate("text", label="ambiguous", y=0.45, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      annotate("text", label="populations", y=0.1, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")
    print(p)
    dev.off()
    
    # plot means and 95% credible interval for GDIs
    dat <- melt(dat)
    dat_mean <- aggregate(dat[,2], list(dat$variable), mean)
    colnames(dat_mean) <- c("species", "mean")
    dat_mean$ci_low <- apply(as.data.frame(allGDIList_byPrior[[i]]), 2, function(x) ci(x, method="HDI", ci=0.95)$CI_low)
    dat_mean$ci_high <- apply(as.data.frame(allGDIList_byPrior[[i]]), 2, function(x) ci(x, method="HDI", ci=0.95)$CI_high)
    colnames(dat_mean) <- c("species", "mean", "ci_low", "ci_high")
		
	print(dat_mean)
        
    pdf(file=paste("priors-tau", prior_df[i,1], prior_df[i,2], "theta", prior_df[i,3], prior_df[i,4], "_means.pdf", sep=""), width=plotWidth, height=plotHeight)
    p <- ggplot(data=dat_mean, aes(x=species, y=mean)) +
      ylim(c(0,1)) +
      geom_hline(yintercept = 0.2, lty=2) +
      geom_hline(yintercept = 0.7, lty=2) +
      geom_errorbar(aes(ymin=ci_low, ymax=ci_high), width=.15) +
      geom_point() +
      labs(y="GDI", x="species")+
      ggtitle(label=paste("priors: tau(", prior_df[i,1], ",", prior_df[i,2], "), theta(", prior_df[i,3], ",", prior_df[i,4], ")", sep="")) +
      annotate("text", label="species", y=0.85, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      annotate("text", label="ambiguous", y=0.45, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      annotate("text", label="populations", y=0.1, x=-Inf, vjust=1.0, hjust=0.5, angle=90) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")
    print(p)
    dev.off()

  }

}
                              
checkConvergence <- function(wd, nreps, burnin) {
  replicate = 0
  setwd(wd)
  
  # get all BPP dirs in the wd
  dirs <- list.dirs(path=wd)
  dirs <- grep("tau", dirs, value=TRUE)
  
  # list to store all gdi estimates
  allMCMCList <- vector(mode = "list", length = length(nreps))
  count=1
  
  # loop through all directories starting with tau
  for(i in dirs) {
    #i = dirs[1]
    replicate = replicate + 1
    rep <- strsplit(i, "/")
    rep <- rep[[1]][length(rep[[1]])]
    
    # read in the posterior
    allMCMCList[[replicate]] <- read.table(paste(i, "/mcmc.txt", sep=""), header=TRUE)[,-1] # drop the first column (number of genes)
	# remove burnin as percentage of the chain
	gen <- nrow(allMCMCList[[replicate]])
	if(burnin > 0){
		burnin_gen <- round(gen*burnin, digits=0)
		allMCMCList[[replicate]] <- allMCMCList[[replicate]][-c(1:burnin_gen),]
	}	
	prior <- strsplit(i, "/")
    # get model name
    model <- prior[[1]][length(prior[[1]])-1]
    # get prior name
    pr <- prior[[1]][length(prior[[1]])]
    pr <- strsplit(pr, "-")
    pr <- Reduce(merge,pr)
    pr <- paste(pr[-length(pr)], collapse="-")
    
    # if we have reached the number of replicates, merge the dataframes in allMCMCList
    if (replicate == nreps){
      message(paste("running diagnostics for ", model, ", priors: ", pr, sep=""))
      replicate = 0
      allmcmc <- ldply(allMCMCList, rbind)
      # convert to mcmc object
      allmcmc <- as.mcmc(allmcmc)
      # convert to mcmc.list object
      allMCMCList <- as.mcmc.list(lapply(allMCMCList, as.mcmc))
      # for each parameter, plot the combined trace
      for(j in 1:ncol(allmcmc)){
        png(file=paste(model, "/", pr, "_", dimnames(allmcmc)[[2]][j], "-trace.png", sep=""), width=8, height=4, units="in", res=300)
        plot(allmcmc[,j])
        dev.off()
      }
      # write effective sample sizes to table
      write.table(effectiveSize(allmcmc), file = paste(model, "/", pr, "_ESS.txt", sep=""), col.names = FALSE)
      # write Galman and Rubin convergence diagnostic (PSRF) to file
      sink(file = paste(model, "/", pr, "_PSRF.txt", sep=""))
      print(gelman.diag(allMCMCList))
      sink()
    }
    count = count + 1
  }
}


