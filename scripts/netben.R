library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
data(dream4)
data(vignette)

# Define a wrapper function
FastBMA <- function(data, threshold=thr){
  edges <- networkBMA(data = data, nTimePoints = nTimePoints,
                      control=fastBMAcontrol(fastgCtrl=fastgControl(optimize=4)))
  g <- graph.data.frame(edges)
  adj <- get.adjacency(g, attr='PostProb',sparse=FALSE)
  
  if (length(colnames(data)) > length(colnames(adj))) {
    missing <- colnames(data)[(!(colnames(data) %in% colnames(adj)))]
    adj <- as.data.frame(adj)
    adj[, missing] <- 0
    adj[missing, ] <- 0
    adj <- as.matrix(adj)
  }
  
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

# Define a wrapper function
ScanBMA <- function(data, threshold=thr){
  edges <- networkBMA(data = data, nTimePoints = nTimePoints)
  
  g <- graph.data.frame(edges)
  adj <- get.adjacency(g, attr='PostProb',sparse=FALSE)

  if (length(colnames(data)) > length(colnames(adj))) {
    missing <- colnames(data)[(!(colnames(data) %in% colnames(adj)))]
    adj <- as.data.frame(adj)
    adj[, missing] <- 0
    adj[missing, ] <- 0
    adj <- as.matrix(adj)
  }
  
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

# Define a wrapper function
BNFinder <- function(data, lim=0, sub=0, k=10, threshold=thr, yeast=yeast.data){
  path <- path_to_bnf
  dat_name <- in.name
  write.table(t(data), paste(path, dat_name, sep=""), sep="\t", quote=FALSE)
  
  if (yeast=="gnw") {
    fConn <- file(paste(path, dat_name, sep=""), 'r+')
    Lines <- readLines(fConn) 
    firstLine <- paste("#regulators", paste(regs, collapse = " "))
    writeLines(c(firstLine, Lines), con = fConn)
    close(fConn)
  } else if (yeast=="bma") {
    fConn <- file(paste(path, dat_name, sep=""), 'r+')
    Lines <- readLines(fConn) 
    u <- unique(referencePairs$Regulator)
    firstLine <- paste("#regulators", paste(u[which(u %in% colnames(data))], collapse = " "))
    writeLines(c(firstLine, Lines), con = fConn)
    close(fConn)
  }

  cores = ""
  if (k!=0) {
    cores = paste(" -k ", k, sep="")
  }
  if ((lim==0)&&(sub==0)) {
    args <- paste(" -e ", path, dat_name, cores, " -v -n ", out.name, sep="")
  } else if ((lim==0)&&(sub!=0)) {
    args <- paste(" -e ", path, dat_name, cores, " -i ", sub, " -v -n ", out.name, sep="")
  } else if ((lim!=0)&&(sub==0)) {
    args <- paste(" -e ", path, dat_name, cores, " -l ",lim, " -v -n ", out.name, sep="")
  } else if ((lim!=0)&&(sub!=0)) {
    args <- paste(" -e ", path, dat_name, cores, " -l ",lim, " -i ", sub, " -v -n ", out.name, sep="")
  }

  cmd <- paste(path, "bnf", args, sep="")
  print(cmd)
  system(cmd)
  
  res <- read.table(out.name)
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

  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
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

aracne.cwrap <- function(data, threshold=thr) {
  adj <- aracne.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

c3net.cwrap <- function(data, threshold=thr) {
  adj <- c3net.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

clr.cwrap <- function(data, threshold=thr) {
  adj <- clr.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

Genie3.cwrap <- function(data, threshold=thr) {
  adj <- Genie3.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

mrnet.cwrap <- function(data, threshold=thr) {
  adj <- mrnet.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

mutrank.cwrap <- function(data, threshold=thr) {
  adj <- mutrank.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

mrnetb.cwrap <- function(data, threshold=thr) {
  adj <- mrnetb.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

pcit.cwrap <- function(data, threshold=thr) {
  adj <- pcit.wrap(data)
  adj <- adj[colnames(data),colnames(data)]
  if (threshold > 0) {
    adj <- adj/max(adj)
    adj[adj < threshold] <- 0
  }
  return(adj)
}

methods.all <- c("ScanBMA", "FastBMA",
                "aracne.cwrap","c3net.cwrap","clr.cwrap",
                "Genie3.cwrap","mrnet.cwrap",
                "mutrank.cwrap","mrnetb.cwrap","pcit.cwrap",
                "BNFinder",
                "BNFinderL3", "BNFinderL2", "BNFinderL1",
                "BNFinderI30", "BNFinderI20", "BNFinderI10", "BNFinderI5",
                "BNFinderL3I30", "BNFinderL3I20", "BNFinderL3I10", "BNFinderL3I5",
                "BNFinderL2I30", "BNFinderL2I20", "BNFinderL2I10", "BNFinderL2I5",
                "BNFinderL1I30", "BNFinderL1I20", "BNFinderL1I10", "BNFinderL1I5")

# Register it to all.fast methods
RegisterWrapper(methods.all)
# Register it to all methods
RegisterWrapper(methods.all, all.fast=FALSE)

benchmark <- function(data, true.net, methods, sym, no.top) {
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

dream.bench <- function(data, gold, methods, sym, no.top) {
  results <- c()
  for (network in 1:1) {
    print(paste("Processing network", network))
    dat <- data[[network]][, -c(1:2)]
    true.net <- gold[[network]]
    true.net <- graph.data.frame(true.net)
    true.net <- get.adjacency(true.net, attr='edge', sparse=FALSE)
    true.net <- as.matrix(true.net[colnames(dat), colnames(dat)])
    result <- benchmark(dat, true.net, methods, sym, no.top)
    results <- c(results, list(result)) 
  }
  
  save(results, file = paste(dname,"top", no.top, sym, ".Rdata", sep=""))
  arr <- abind(results, along=3)
  average <- rowMeans(arr, dims = 2)
  average <- as.data.frame(t(average))
  Time <- rowMeans(average[,c("TimeRep1", "TimeRep2", "TimeRep3")])
  average <- cbind(average, Time)
  average <- average[,-which(colnames(average) %in% c("TimeRep1", "TimeRep2", "TimeRep3"))]  
  return(average)    
}
 
yeast.bench <- function(data, gold, methods, sym, no.top) {
  true.net <- gold
  true.net$edge <- 1
  true.net <- graph.data.frame(true.net)
  true.net <- get.adjacency(true.net, attr='edge', sparse=FALSE)
  true.net <- as.matrix(true.net[colnames(data), colnames(data)])
  result <- benchmark(data, true.net, methods, sym, no.top)

  result <- as.data.frame(t(result))
  Time <- rowMeans(result[,c("TimeRep1", "TimeRep2", "TimeRep3")])
  result <- cbind(result, Time)
  return(result)    
}


path_to_bnf = "/home/aln/Science/URGENT/article/article-bnf-suppl/bnfinder/"

## Dream4 tests

methods <- methods.all

#This is a shame, but I use global variables:)
data = dream4ts10
dname = "dream4ts10"
in.name = "dream10.in"
out.name = "dream10.out"
gold = dream4gold10
nTimePoints <- length(unique(data[[1]]$time))
yeast.data="none"

thr=0
result <- dream.bench(data, gold, methods, no.top=20, sym=FALSE)
write.table(result, paste(dname,"top20.txt", sep=""), sep="\t", quote=FALSE)

result <- dream.bench(data, gold, methods, no.top=20, sym=TRUE)
write.table(result, paste(dname,"top20sym.txt", sep=""), sep="\t", quote=FALSE)

thr=0.5
result <- dream.bench(data, gold, methods, no.top=100, sym=FALSE)
write.table(result, paste(dname,"top100.txt", sep=""), sep="\t", quote=FALSE)

result <- dream.bench(data, gold, methods, no.top=100, sym=TRUE)
write.table(result, paste(dname,"top100sym.txt", sep=""), sep="\t", quote=FALSE)

#Evaluate


ngenes <- dim(syntren300.data)[2]
nlinks <- ngenes^2-ngenes
no.topedges=100
no.edges <- round(nlinks*no.topedges/100)

net <- clr.wrap(syntren300.data)

r <- evaluate(net,syntren300.net,extend=no.edges,sym=FALSE)
aupr(r)
r <- evaluate(net,syntren300.net,extend=0,sym=FALSE)
aupr(r)

net <- net/max(net)
net[net<0.5] <- 0

syntren300.net.nozeros <- syntren300.net
syntren300.net.nozeros[syntren300.net.nozeros==0] <- 0.9
syntren300.net.nozeros[1,1] <- 0
r <- evaluate(syntren300.net.nozeros,syntren300.net,extend=no.edges,sym=FALSE)
aupr(r)
auroc(r)

r <- evaluate(syntren300.net.nozeros,syntren300.net,extend=0,sym=FALSE)
aupr(r)
auroc(r)

t <- minet::validate(syntren300.net.nozeros,syntren300.net)
minet::auc.pr(t)
minet::auc.roc(t)

data(dream4)
network <- 2
reference <- dream4gold10[[network]]
nGenes <- length(unique(c(reference[,1],reference[,2])))
nPossibleEdges <- nGenes^2
reference <- reference[reference[,3] == 1,1:2]
nTimePoints <- length(unique(dream4ts10[[network]]$time))
edges1ts10 <- networkBMA( data = dream4ts10[[network]][,-(1:2)],
                          nTimePoints = nTimePoints, prior.prob = 0.1,
                          self = FALSE)
size <- nPossibleEdges - nGenes
thr=0
yeast.data=""
bnf <- BNFinderI10(dream4ts10[[network]][,-(1:2)])
g <- graph.adjacency(bnf, weighted=TRUE, mode="directed")
edges <- as.data.frame(get.edgelist(g))
edges$weight <- E(g)$weight


contingencyTables <- contabs(network = edges1ts10, reference = reference,
                             size = size, thresholds = 0.1)
scores(contingencyTables, what = c("sensitivity", "specificity", "FDR"))
roc(contingencyTables, plotit = TRUE)
prc(contingencyTables, plotit = TRUE)

colnames(edges) <- names(edges1ts10)
contingencyTables <- contabs(network = edges, reference = reference,
                             size = size, thresholds = 0.1)
scores(contingencyTables, what = c("sensitivity", "specificity", "FDR"))
roc(contingencyTables, plotit = TRUE)
prc(contingencyTables, plotit = TRUE)


## Brem data tests

yeast.data="bma"
data = as.data.frame(t(brem.data))
gold = referencePairs
dname = "brem"
in.name = "brem.in"
out.name = "brem.out"

methods <- c("aracne.cwrap","c3net.cwrap","clr.cwrap",
                            "Genie3.cwrap","mrnet.cwrap",
                            "mutrank.cwrap","mrnetb.cwrap","pcit.cwrap",
                            "BNFinderL1I10")
                            
result <- yeast.bench(data, gold, methods, no.top=20, sym=FALSE)

# Priors for time series
g <- graph.adjacency(1-reg.prob, weighted=TRUE, mode="directed")
edges <- as.data.frame(get.edgelist(g))
edges$weight <- E(g)$weight
edges <- edges[,c(2, 3, 1)]
edges <- edges[which(edges$weight<0.9),]
known <- reg.known[, c(1, 2)]
known$prob <- 0.2
known <- known[, c(2, 3, 1)]
colnames(edges) <- colnames(known)
edges <- rbind(edges, known)
edges$pe <- "#prioredge"
edges <- edges[,c(4, 1, 2, 3)]
c <- apply(edges, 1, function(x) {paste(x, collapse=" ")})



true.net <- dream4gold10[[1]]
true.net <- graph.data.frame(true.net)
true.net <- get.adjacency(true.net, attr='edge', sparse=FALSE)
true.net <- as.matrix(true.net[colnames(dream4ts10[[1]][,-c(1:2)]), colnames(dream4ts10[[1]][,-c(1:2)])])
data <- dream4ts10[[1]][,-c(1:2)]
bnf <- BNFinderL3I5(data)
bma <- FastBMA(data)
clr <- clr.wrap(data)

aupr <- netbenchmark.data(data=dream4ts10[[1]][,-c(1:2)],
                          methods=c("BNFinderL3I5","FastBMA"), sym = FALSE,
                          no.topedges = 50, true.net = true.net, plot=TRUE)



library(minet)
syntren300.net.nozeros <- syntren300.net
length(which(syntren300.net.nozeros==0))
syntren300.net.nozeros[syntren300.net.nozeros==0] <- runif(89532, 0.0, 1.1)

bf <- minet::validate(bnf,true.net)
bm <- minet::validate(bma,true.net)
cl <- minet::validate(clr, true.net)

auc.pr(bf)
auc.pr(bm)
auc.pr(cl)
auc.roc(bf)
auc.roc(bm)
auc.roc(cl)

dev <- minet::show.pr(bf, col="green", type="b")
dev <- minet::show.pr(bm, col="blue", type="b", dev = dev)
minet::show.pr(cl, device=dev, col="red",type="b")

dev.off()
dev <- minet::show.roc(bf, col="green", type="b")
dev <- minet::show.roc(bm, col="blue", type="b", dev = dev)
minet::show.roc(cl, device=dev, col="red",type="b")
