library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
data(dream4)
data(vignette)


# Read Yeast TFs
#Yeast.TFs <- read.table("../data/RegulationTwoColumnTable_Documented_2013927.tsv",
#                        check.names = FALSE, header=FALSE, sep=";", quote="")
#Yeast.TFs <- levels(unique(Yeast.TFs[,1]))
#gene_to_ORF <- select(org.Sc.sgd.db, keys=Yeast.TFs, columns = c("ORF"),
#                      keytype="COMMON")


# Define a wrapper function
ScanBMA <- function(data){
  edges <- networkBMA(data = data, nTimePoints = nrow(data))
  g <- graph.data.frame(edges)
  adj <- get.adjacency(g, attr='PostProb',sparse=FALSE)
  return(adj[colnames(data),colnames(data)])
}

# Define a wrapper function
BNFinder <- function(data, lim=0, sub=0, k=10){
  path <- ""
  dat_name <- "input.txt"
  write.table(t(data), paste(path, dat_name, sep=""), sep="\t", quote=FALSE)
  cores = ""
  if (k!=0) {
    cores = paste(" -k ", k, sep="")
  }
  if ((lim==0)&&(sub==0)) {
    args <- paste(" -e ", path, dat_name, cores, " -v -t output.txt -n output.tsv", sep="")
  } else if ((lim==0)&&(sub!=0)) {
    args <- paste(" -e ", path, dat_name, cores, " -i ", sub, " -v -t output.txt -n output.tsv", sep="")
  } else if ((lim!=0)&&(sub==0)) {
    args <- paste(" -e ", path, dat_name, cores, " -l ",lim, " -v -t output.txt -n output.tsv", sep="")
  } else if ((lim!=0)&&(sub!=0)) {
    args <- paste(" -e ", path, dat_name, cores, " -l ",lim, " -i ", sub, " -v -t output.txt -n output.tsv", sep="")
  }

  cmd <- paste(path, "bnf", args, sep="")
  print(cmd)
  system(cmd)
  
  res <- read.table("output.tsv")
  res <- res[, c(1, 3, 2)]
  if (sub==0) {
    res[,3] <- 1
  }

  g <- graph.data.frame(res)
  adj <- get.adjacency(g, attr='V2',sparse=FALSE)
  if (length(colnames(data)) > length(colnames(adj))) {
    missing <- colnames(data)[(!(colnames(data) %in% colnames(adj)))]
    adj <- as.data.frame(adj)
    adj[, missing] <- 0
    adj[missing, ] <- 0
    adj <- as.matrix(adj)
  }
  
  return(adj[colnames(data),colnames(data)])
}

BNFinderL3 <- function(data){
  return(BNFinder(data, lim=3, sub=0))
}

BNFinderL2 <- function(data){
  return(BNFinder(data, lim=2, sub=0))
}

BNFinderL1 <- function(data){
  return(BNFinder(data, lim=1, sub=0))
}

BNFinderI5 <- function(data){
  return(BNFinder(data, lim=0, sub=5))
}

BNFinderI10 <- function(data){
  return(BNFinder(data, lim=0, sub=10))
}

BNFinderI20 <- function(data){
  return(BNFinder(data, lim=0, sub=20))
}

BNFinderI30 <- function(data){
  return(BNFinder(data, lim=0, sub=30))
}

BNFinderL3I30 <- function(data){
  return(BNFinder(data, lim=3, sub=30))
}

BNFinderL3I20 <- function(data){
  return(BNFinder(data, lim=3, sub=20))
}

BNFinderL3I10 <- function(data){
  return(BNFinder(data, lim=3, sub=10))
}

BNFinderL3I5 <- function(data){
  return(BNFinder(data, lim=3, sub=5))
}

BNFinderL2I30 <- function(data){
  return(BNFinder(data, lim=2, sub=30))
}

BNFinderL2I20 <- function(data){
  return(BNFinder(data, lim=2, sub=20))
}

BNFinderL2I10 <- function(data){
  return(BNFinder(data, lim=2, sub=10))
}

BNFinderL2I5 <- function(data){
  return(BNFinder(data, lim=2, sub=5))
}

BNFinderL1I30 <- function(data){
  return(BNFinder(data, lim=1, sub=30))
}

BNFinderL1I20 <- function(data){
  return(BNFinder(data, lim=1, sub=20))
}

BNFinderL1I10 <- function(data){
  return(BNFinder(data, lim=1, sub=10))
}

BNFinderL1I5 <- function(data){
  return(BNFinder(data, lim=1, sub=5))
}


# Register it to all.fast methods
RegisterWrapper(c("ScanBMA", "BNFinder",
                  "BNFinderL3", "BNFinderL2", "BNFinderL1",
                  "BNFinderI30", "BNFinderI20", "BNFinderI10", "BNFinderI5",
                  "BNFinderL3I30", "BNFinderL3I20", "BNFinderL3I10", "BNFinderL3I5",
                  "BNFinderL2I30", "BNFinderL2I20", "BNFinderL2I10", "BNFinderL2I5",
                  "BNFinderL1I30", "BNFinderL1I20", "BNFinderL1I10", "BNFinderL1I5"))
                  # Register it to all methods
RegisterWrapper(c("ScanBMA", "BNFinder",
                  "BNFinderL3", "BNFinderL2", "BNFinderL1",
                  "BNFinderI30", "BNFinderI20", "BNFinderI10", "BNFinderI5",
                  "BNFinderL3I30", "BNFinderL3I20", "BNFinderL3I10", "BNFinderL3I5",
                  "BNFinderL2I30", "BNFinderL2I20", "BNFinderL2I10", "BNFinderL2I5",
                  "BNFinderL1I30", "BNFinderL1I20", "BNFinderL1I10", "BNFinderL1I5"),
                all.fast=FALSE)

benchmark <- function(data, true.net, methods, sym=TRUE, no.top=50) {
  auroc <- netbenchmark.data(data=data, eval="AUROC", methods=methods, sym = sym,
                             no.topedges = no.top, true.net = true.net)
  aupr <- netbenchmark.data(data=data, methods=methods, sym = sym,
                          no.topedges = no.top, true.net = true.net)
  tp <- netbenchmark.data(data=data, eval="no.truepos", methods=methods, sym = sym,
                           no.topedges = no.top, true.net = true.net)
      
  result <- rbind(auroc[[1]], aupr[[1]], tp[[1]])
  cpu <- rbind(auroc[[2]], aupr[[2]], tp[[2]])
  #cpu <- as.data.frame(t(colMeans(cpu)))
  cpu$rand <- 0
  result <- rbind(result, cpu)
  rownames(result) <- c("AUROC", "AUPR", "TP", "TimeRep1", "TimeRep2", "TimeRep3")
  return(result)
}

methods <- c("ScanBMA", "BNFinder", "aracne.wrap","c3net.wrap","clr.wrap",
                 "Genie3.wrap","mrnet.wrap",
                 "mutrank.wrap","mrnetb.wrap","pcit.wrap")
    
methods <- c("ScanBMA", "aracne.wrap","c3net.wrap","clr.wrap",
             "Genie3.wrap","mrnet.wrap",
             "mutrank.wrap","mrnetb.wrap","pcit.wrap",
	     "BNFinder",
             "BNFinderL3", "BNFinderL2", "BNFinderL1",
             "BNFinderI30", "BNFinderI20", "BNFinderI10", "BNFinderI5",
	     "BNFinderL3I30", "BNFinderL3I20", "BNFinderL3I10", "BNFinderL3I5",
             "BNFinderL2I30", "BNFinderL2I20", "BNFinderL2I10", "BNFinderL2I5",
             "BNFinderL1I30", "BNFinderL1I20", "BNFinderL1I10", "BNFinderL1I5")

methods <- c("ScanBMA", "BNFinder")
#methods <- c("GeneNet.wrap")

dream.bench <- function(data, gold, methods, sym=TRUE, no.top=50) {
  results <- c()
  for (network in 1:length(data)) {
    print(paste("Processing network", network))
    dat <- data[[network]][, -c(1:2)]
    true.net <- gold[[network]]
    true.net <- graph.data.frame(true.net)
    true.net <- get.adjacency(true.net, attr='edge',sparse=FALSE)
    true.net <- as.matrix(true.net[colnames(dat), colnames(dat)])
    result <- benchmark(dat, true.net, methods, sym, no.top)
    results <- c(results, list(result)) 
  }
  
  save(results, file = paste("top", no.top, sym, ".Rdata", sep=""))
  arr <- abind(results, along=3)
  average <- rowMeans(arr, dims = 2)
  average <- as.data.frame(t(average))
  Time <- rowMeans(average[,c("TimeRep1", "TimeRep2", "TimeRep3")])
  average <- cbind(average, Time)
  average <- average[,-which(colnames(average) %in% c("TimeRep1", "TimeRep2", "TimeRep3"))]  
  return(average)    
}
 
data = dream4ts10
gold = dream4gold10

result <- dream.bench(data, gold, methods, no.top=20, sym=FALSE)
write.table(result, "top20.txt", sep="\t", quote=FALSE)
result <- dream.bench(data, gold, methods, no.top=50, sym=FALSE)
write.table(result, "top50.txt", sep="\t", quote=FALSE)
result <- dream.bench(data, gold, methods, no.top=80, sym=FALSE)
write.table(result, "top80.txt", sep="\t", quote=FALSE)
result <- dream.bench(data, gold, methods, no.top=100, sym=FALSE)
write.table(result, "top100.txt", sep="\t", quote=FALSE)


result <- dream.bench(data, gold, methods, no.top=20, sym=TRUE)
write.table(result, "top20sym.txt", sep="\t", quote=FALSE)
result <- dream.bench(data, gold, methods, no.top=50, sym=TRUE)
write.table(result, "top50sym.txt", sep="\t", quote=FALSE)
result <- dream.bench(data, gold, methods, no.top=80, sym=TRUE)
write.table(result, "top80sym.txt", sep="\t", quote=FALSE)
result <- dream.bench(data, gold, methods, no.top=100, sym=TRUE)
write.table(result, "top100sym.txt", sep="\t", quote=FALSE)


#View(result)

