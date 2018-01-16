library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
library(minet)
data(dream4)
data(vignette)


BNFinder <- function(data, data.name, regulators=NULL, priors=NULL, lim=0, sub=0, k=0, path=""){
  in.name <- paste(getwd(), "/../bnf_in/", data.name, ".in", sep="")
  write.table(t(data), in.name, sep="\t", quote=FALSE)

  if (!is.null(priors)) {
    priors$pe <- "#prioredge"
    priors <- priors[,c(4, 2, 3, 1)]
    p <- apply(priors, 1, function(x) {paste(x, collapse=" ")})
    fConn <- file(in.name, 'r+')
    Lines <- readLines(fConn) 
    writeLines(c(p, Lines), con = fConn) 
    close(fConn)
  }
  
  if (!is.null(regulators)) {
    fConn <- file(in.name, 'r+')
    Lines <- readLines(fConn) 
    firstLine <- paste("#regulators", paste(regulators, collapse = " "))
    writeLines(c(firstLine, Lines), con = fConn) 
    close(fConn)
  }

  cores = ""
  if (k!=0) {
    cores = paste(" -k ", k, sep="")
  }
  
  out.name = ""
  if ((lim==0)&&(sub==0)) {
    out.name = paste(getwd(), "/../bnf_out/", data.name, ".out", sep="")
    args <- paste(" -e ", in.name, cores, " -v -n ", out.name, sep="")
  } else if ((lim==0)&&(sub!=0)) {
    out.name = paste(getwd(), "/../bnf_out/", data.name, "I", sub, ".out", sep="")
    args <- paste(" -e ", in.name, cores, " -i ", sub, " -v -n ", out.name, sep="")
  } else if ((lim!=0)&&(sub==0)) {
    out.name = paste(getwd(), "/../bnf_out/", data.name, "L", lim, ".out", sep="")
    args <- paste(" -e ", in.name, cores, " -l ",lim, " -v -n ", out.name, sep="")
  } else if ((lim!=0)&&(sub!=0)) {
    out.name = paste(getwd(), "/../bnf_out/", data.name, "L", lim, "I", sub, ".out", sep="")
    args <- paste(" -e ", in.name, cores, " -l ",lim, " -i ", sub, " -v -n ", out.name, sep="")
  }
  
  cmd <- paste(path, "bnf", args, sep="")
  print(cmd)
  system(cmd)
  
  res <- read.table(out.name)
  res <- res[, c(1, 3, 2)]
  if (sub==0) {
    res[,3] <- 1
  }
  colnames(res) <- c("Regulator", "Target", "Weight")
  if (k>0) {
    time <- read.table(paste(getwd(), "/time", k, ".txt", sep=""))
  }
  return(list(time, res))
}

params <- data.frame(lim=c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3),
                     sub = c(0, 5, 10, 20, 30, 40, 0, 5, 10, 20, 30, 40, 0, 5, 10, 20, 30, 40, 0, 5, 10, 20, 30, 40))

path = paste(getwd(), "/../bnfinder/", sep="")

### Dream4 10 and 100 TS

DREAM = FALSE

if (DREAM) {

  data = dream4ts10
  gold = dream4gold10
  data.name = "dream4ts10"

  results <- c()
  nn <- length(data)
  for (network in 1:nn) {
    print(paste("Processing network", network))
    dat <- data[[network]][, -c(1:2)]
    for (i in 1:nrow(params)) {
      net.name <- paste(data.name, "n", network, sep="")
      result <- BNFinder(dat, net.name, lim=params[i, 1], sub=params[i, 2], k=2, path=path)
      results <- c(results,
                   list(list(network=net.name,params=paste("L",params[i,1],"I",params[i,2],sep=""),
                             time.sec=as.numeric(result[[1]]), inf.net=result[[2]])))
    }
  }  

}


### Yeast

# Read Yeast TFs
Yeast.TFs <- read.table("../data/RegulationTwoColumnTable_Documented_2013927.tsv",
                        check.names = FALSE, header=FALSE, sep=";", quote="")
Yeast.TFs <- levels(unique(Yeast.TFs[,1]))
gene_to_ORF <- select(org.Sc.sgd.db, keys=Yeast.TFs, columns = c("ORF"),
                      keytype="COMMON")

# Priors
g <- graph.adjacency(1-reg.prob, weighted=TRUE, mode="directed")
edges <- as.data.frame(get.edgelist(g))
edges$weight <- E(g)$weight
edges <- edges[which(edges$weight<0.9),]
known <- reg.known[, c(1, 2)]
known$prob <- 0.2
colnames(edges) <- colnames(known)
edges <- rbind(edges, known)


## Brem

BREM=FALSE
if (BREM) {
  data = t(brem.data)
  gold = referencePairs
  data.name = "brem"
  
  regs <- colnames(data)[which(colnames(data) %in% gene_to_ORF$ORF)]
  results <- c()
  
  params <- data.frame(lim=c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
                       sub = c(0, 5, 10, 20, 30, 0, 5, 10, 20, 30, 0, 5, 10, 20, 30))
  
  for (i in 1:nrow(params)) {
    result <- BNFinder(data, data.name, regulators=regs, priors=edges,
                      lim=params[i,1], sub=params[i,2], k=1, path=path)
  
    results <- c(results, list(list(network=data.name,params=paste("L",params[i,1],"I",params[i,2],sep=""),
                                    time.sec=as.numeric(result[[1]]), inf.net=result[[2]])))
  }

}

## Time Series Yeast

data = timeSeries[, -c(1:2)]
gold = referencePairs
data.name = "YeastTS"

results <- c()

params <- data.frame(lim=c(1),
                     sub = c(2))

for (i in 1:nrow(params)) {
  result <- BNFinder(data, data.name, priors=edges,
                     lim=params[i,1], sub=params[i,2], k=1, path=path)
  
  results <- c(results, list(list(network=data.name,params=paste("L",params[i,1],"I",params[i,2],sep=""),
                                  time.sec=as.numeric(result[[1]]), inf.net=result[[2]])))
}



########
save(results, file=paste(data.name, "_results.RData", sep=""))

