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

FastBMA <- function(data, nTimePoints){
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
  adj <- adj/max(adj)
  return(adj)
}

data = dream4ts10
gold = dream4gold10
data.name = "dream4ts10"

load(paste("../bnf_out/", data.name, "_results.RData", sep=""))

eval20top.sym <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))
eval20top <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))

for (i in 1:length(results)) {
  network = as.numeric(substr(results[[i]]$network, nchar(results[[i]]$network), nchar(results[[i]]$network)))
  results[[i]]$inf.net$Weight <- as.numeric(results[[i]]$inf.net$Weight)

  method <- paste("BNFinder", results[[i]]$params, sep="")
  nGenes <- ncol(data[[network]][, -c(1:2)])
  nEdges = round((nGenes^2 - nGenes)*0.2)
  #nEdges = sum(gold[[network]]$edge)
  
  e <- evaluate(EdgeToAdj(results[[i]]$inf.net, colnames(data[[network]][, -c(1:2)]), attr="Weight"),
                EdgeToAdj(gold[[network]], colnames(data[[network]][, -c(1:2)]), attr="edge"),
                sym=FALSE, extend=nEdges)
  eval20top <- rbind(eval20top, c(method, auroc(e, nEdges), aupr(e, nEdges), network, results[[i]]$time.sec))

  e <- evaluate(EdgeToAdj(results[[i]]$inf.net, colnames(data[[network]][, -c(1:2)]), attr="Weight"),
                EdgeToAdj(gold[[network]], colnames(data[[network]][, -c(1:2)]), attr="edge"),
                sym=TRUE, extend=nEdges)
  eval20top.sym <- rbind(eval20top.sym, c(method, auroc(e, nEdges), aupr(e, nEdges), network, results[[i]]$time.sec))
}

eval20top <- eval20top[-1, ]
eval20top <- eval20top[grep("L0I10", eval20top$method),]
eval20top.sym <- eval20top.sym[-1, ]
eval20top.sym <- eval20top.sym[grep("L0I10", eval20top.sym$method),]

methods <- c("FastBMA", "aracne","c3net","clr",
             "Genie3","mrnet",
             "mutrank","mrnetb","pcit")
functions <- c("FastBMA", "aracne.wrap","c3net.wrap","clr.wrap",
               "Genie3.wrap","mrnet.wrap",
               "mutrank.wrap","mrnetb.wrap","pcit.wrap")

methods.nonsym <- c("FastBMA", "Genie3")
functions.nonsym <- c("FastBMA", "Genie3.wrap")

for (i in 1:length(data)) {
  network=i
  d <- data[[network]][, -c(1,2)]
  nGenes <- ncol(data[[network]][, -c(1:2)])
  nEdges = round((nGenes^2 - nGenes)*0.2)
  #nEdges = sum(gold[[network]]$edge)
  nTimePoints <- length(unique(data[[network]]$time))
  
  for (j in 1:length(methods)) {
    method = methods[j]
    start_time <- Sys.time()
    if (method=="FastBMA") {
      inf <- do.call(functions[j], list(d, nTimePoints))
    } else {
      inf <- do.call(functions[j], list(d))
    }
    end_time <- Sys.time()
    time <- as.numeric(end_time - start_time)
  
    e <- evaluate(inf, EdgeToAdj(gold[[network]], colnames(data[[network]][, -c(1:2)]), attr="edge"),
                  sym=TRUE, extend=nEdges)
    eval20top.sym <- rbind(eval20top.sym, c(method, auroc(e, nEdges), aupr(e, nEdges), network, time))
  }
  
  for (j in 1:length(methods.nonsym)) {
    method = methods.nonsym[j]
    start_time <- Sys.time()
    if (method=="FastBMA") {
      inf <- do.call(functions[j], list(d, nTimePoints))
    } else {
      inf <- do.call(functions[j], list(d))
    }
    end_time <- Sys.time()
    time <- as.numeric(end_time - start_time)

    e <- evaluate(inf, EdgeToAdj(gold[[network]], colnames(data[[network]][, -c(1:2)]), attr="edge"),
                  sym=FALSE, extend=nEdges)
    print(e)    
    eval20top <- rbind(eval20top, c(method, auroc(e, nEdges), aupr(e, nEdges), network, time))
  }
  
}

eval20top <- eval20top[order(eval20top$no.net),]
eval20top.sym <- eval20top.sym[order(eval20top.sym$no.net),]
