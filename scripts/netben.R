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
  return(adj[,colnames(data)])
}

# Define a wrapper function
BNFinder <- function(data){
  path <- ""
  dat_name <- "input.txt"
  write.table(t(data), paste(path, dat_name, sep=""), sep="\t", quote=FALSE)
  args <- paste(" -e ", path, dat_name, " -v -t output.txt -n output.tsv", sep="")
  system(paste(path, "bnf", args, sep=""))
  
  res <- read.table("output.tsv")
  res <- res[, c(1, 3, 2)]
  res[,3] <- 1
  g <- graph.data.frame(res)
  adj <- get.adjacency(g, attr='V2',sparse=FALSE)
  if (length(colnames(data)) > length(colnames(adj))) {
    missing <- colnames(data)[(!(colnames(data) %in% colnames(adj)))]
    adj <- as.data.frame(adj)
    adj[, missing] <- 0
    adj[missing, ] <- 0
    adj <- as.matrix(adj)
  }
  
  return(adj[,colnames(data)])
}

# Register it to all.fast methods
RegisterWrapper(c("ScanBMA", "BNFinder"))
# Register it to all methods
RegisterWrapper(c("ScanBMA", "BNFinder"), all.fast=FALSE)

benchmark <- function(data, true.net, methods, sym=TRUE, no.top=50) {
  auroc <- netbenchmark.data(data=data, eval="AUROC", methods=methods, sym = sym,
                             no.topedges = no.top, true.net = true.net)
  aupr <- netbenchmark.data(data=data, methods=methods, sym = sym,
                          no.topedges = no.top, true.net = true.net)
  tp <- netbenchmark.data(data=data, eval="no.truepos", methods=methods, sym = sym,
                           no.topedges = no.top, true.net = true.net)
      
  result <- rbind(auroc[[1]], aupr[[1]], tp[[1]])
  cpu <- rbind(auroc[[2]], aupr[[2]], tp[[2]])
  cpu <- as.data.frame(t(colMeans(cpu)))
  cpu$rand = 0
  result <- rbind(result, cpu)
  rownames(result) <- c("AUROC", "AUPR", "TP", "Time")
  return(result)
}

methods <- c("ScanBMA", "BNFinder", "aracne.wrap","c3net.wrap","clr.wrap",
                 "Genie3.wrap","mrnet.wrap",
                 "mutrank.wrap","mrnetb.wrap","pcit.wrap")
    
#methods <- c("BNFinder")
#methods <- c("clr.wrap", "ScanBMA")
#methods <- c("GeneNet.wrap")

dream.bench <- function(data, gold, methods, sym=TRUE) {
  results <- c()
  for (network in 2:2) {
    dat <- data[[network]][, -c(1:2)]
    true.net <- gold[[network]]
    true.net <- graph.data.frame(true.net)
    true.net <- get.adjacency(true.net, attr='edge',sparse=FALSE)
    true.net <- as.matrix(true.net[colnames(dat), colnames(dat)])
    result <- benchmark(dat, true.net, methods, sym)
    results <- c(results, list(result)) 
  }
  
  arr <- abind(results, along=3)
  average <- rowMeans(arr, dims = 2)
  
  return(average)    
}
  
result <- dream.bench(dream4ts10, dream4gold10, methods, sym=TRUE)
write.table(result, "result.txt", sep="\t", quote=FALSE)

#View(result)

