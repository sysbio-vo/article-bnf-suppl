library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
library(minet)
data(dream4)
data(vignette)

EdgeToAdj <- function(edgeList, colnames, attr) {
  g <- graph.data.frame(edgeList)
  adj <- get.adjacency(g, attr=attr,sparse=FALSE)
  if (length(colnames) > length(colnames(adj))) {
    missing <- colnames[(!(colnames %in% colnames(adj)))]
    adj <- as.data.frame(adj)
    adj[, missing] <- 0
    adj[missing, ] <- 0
    adj <- as.matrix(adj)
  }
  
  adj <- adj[colnames,colnames]
  return(adj)
}

FastBMA <- function(data, nTimePoints, priors = NULL, known = NULL){
  edges <- networkBMA(data = data, nTimePoints = nTimePoints, prior.prob = priors, known = known,
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
  adj <- adj/max(adj)
  return(adj)
}

data = t(brem.data)
gold = referencePairs
gold$edge <- 1
data.name = "brem"

load(paste("../bnf_out/", data.name, "_results.RData", sep=""))

eval20top.sym <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))
eval20top <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))

for (i in 1:length(results)) {
  results[[i]]$inf.net$Weight <- as.numeric(results[[i]]$inf.net$Weight)
  
  method <- paste("BNFinder", results[[i]]$params, sep="")
  nGenes <- ncol(data)
  nEdges = round((nGenes^2 - nGenes)*0.2)
  #nEdges = sum(gold[[network]]$edge)
  
  e <- evaluate(EdgeToAdj(results[[i]]$inf.net, colnames(data), attr="Weight"),
                EdgeToAdj(gold, colnames(data), attr="edge"), sym=FALSE, extend=nEdges)
  eval20top <- rbind(eval20top, c(method, auroc(e, nEdges), aupr(e, nEdges), data.name, results[[i]]$time.sec))
  
  e <- evaluate(EdgeToAdj(results[[i]]$inf.net, colnames(data), attr="Weight"),
                EdgeToAdj(gold, colnames(data), attr="edge"), sym=TRUE, extend=nEdges)
  eval20top.sym <- rbind(eval20top.sym, c(method, auroc(e, nEdges), aupr(e, nEdges), data.name, results[[i]]$time.sec))
}

eval20top <- eval20top[-1, ]
eval20top.sym <- eval20top.sym[-1, ]

eval20top <- eval20top[grep("L0I10", eval20top$method),]
eval20top.sym <- eval20top.sym[grep("L0I10", eval20top.sym$method),]

methods <- c("aracne","c3net","clr",
             "Genie3","mrnet",
             "mutrank","mrnetb","pcit")
functions <- c("aracne.wrap","c3net.wrap","clr.wrap",
               "Genie3.wrap","mrnet.wrap",
               "mutrank.wrap","mrnetb.wrap","pcit.wrap")

methods.nonsym <- c("Genie3")
functions.nonsym <- c("Genie3.wrap")

nGenes <- ncol(data)
nEdges = round((nGenes^2 - nGenes)*0.2)

for (j in 1:length(methods)) {
  method = methods[j]
  start_time <- Sys.time()
  if (method=="FastBMA") {
    inf <- do.call(functions[j], list(data, nTimePoints, priors = reg.prob, known = reg.known))
  } else {
    inf <- do.call(functions[j], list(d))
  }
  end_time <- Sys.time()
  time <- as.numeric(end_time - start_time)
  
  e <- evaluate(inf, EdgeToAdj(gold, colnames(data), attr="edge"), sym=TRUE, extend=nEdges)
  eval20top.sym <- rbind(eval20top.sym, c(method, auroc(e, nEdges), aupr(e, nEdges), data.name, time))
}

for (j in 1:length(methods.nonsym)) {
  method = methods.nonsym[j]
  start_time <- Sys.time()
  if (method=="FastBMA") {
    inf <- do.call(functions[j], list(data, nTimePoints, priors = reg.prob, known = reg.known))
  } else {
    inf <- do.call(functions[j], list(d))
  }
  end_time <- Sys.time()
  time <- as.numeric(end_time - start_time)
  
  e <- evaluate(inf, EdgeToAdj(gold, colnames(data), attr="edge"), sym=FALSE, extend=nEdges)
  eval20top <- rbind(eval20top, c(method, auroc(e, nEdges), aupr(e, nEdges), data.name, time))
}

