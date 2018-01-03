library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
data(dream4)
data(vignette)


# Read Yeast TFs
Yeast.TFs <- read.table("../data/RegulationTwoColumnTable_Documented_2013927.tsv",
                        check.names = FALSE, header=FALSE, sep=";", quote="")
Yeast.TFs <- levels(unique(Yeast.TFs[,1]))
gene_to_ORF <- select(org.Sc.sgd.db, keys=Yeast.TFs, columns = c("ORF"),
                      keytype="COMMON")


# Define a wrapper function
ScanBMA <- function(data){
  edges <- networkBMA(data = data, nTimePoints = nrow(data))
  g <- graph.data.frame(edges)
  adj <- get.adjacency(g, attr='PostProb',sparse=FALSE)
  return(adj[,colnames(data)])
}

# Register it to all.fast methods
RegisterWrapper("ScanBMA")
# Register it to all methods
RegisterWrapper("ScanBMA", all.fast=FALSE)

benchmark <- function(data, true.net, methods, sym=TRUE, no.top=100) {
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

methods <- c("ScanBMA", "aracne.wrap","c3net.wrap","clr.wrap",
                 "Genie3.wrap","mrnet.wrap",
                 "mutrank.wrap","mrnetb.wrap","pcit.wrap")
    
#methods <- c("clr.wrap", "ScanBMA")
#methods <- c("GeneNet.wrap")

dream.bench <- function(data, gold, methods) {
  results <- c()
  for (network in 1:length(data)) {
    dat <- data[[network]][, -c(1:2)]
    true.net <- gold[[network]]
    true.net <- graph.data.frame(true.net)
    true.net <- get.adjacency(true.net, attr='edge',sparse=FALSE)
    true.net <- as.matrix(true.net[colnames(dat), colnames(dat)])
    result <- benchmark(dat, true.net, methods)
    results <- c(results, list(result)) 
  }

  return(results)
  arr <- abind(results, along=3)
  average <- rowMeans(arr, dims = 2)
  
  return(average[[1]])    
}
  
result <- dream.bench(dream4ts10, dream4gold10, methods)
View(result)

system("./test.sh")
