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

eval_minet <- function(inf, gold, sym) {
  if (sym==FALSE) {
    e <- validate(inf, gold)
  } else {
    if (!isSymmetrix(gold)) {
      gold <- pmax(gold,t(gold))
    }
    e <- validate(inf, gold)
  }
  roc <- auc.roc(e)
  pr <- auc.pr(e)
  
  return(c(roc, pr))
}

eval_netb <- function(inf, gold, sym, no.top) {
  nEdges <- round((ncol(gold)^2 - ncol(gold))*no.top)
  if (sym==FALSE) {
    e <- evaluate(inf, gold, sym=FALSE, extend=nEdges)
  } else {
    e <- evaluate(inf, gold, sym=TRUE, extend=nEdges)
  }
  roc <- auroc(e, nEdges)
  pr <- aupr(e, nEdges)
  
  return(c(roc, pr))
}

eval_methods <- c("netb", "mnet")
emethod <- eval_methods[2]

data = dream4ts100
gold = dream4gold100
data.name = "dream4ts100"

load(paste("../bnf_out/", data.name, "_results.RData", sep=""))

eval.sym <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))
eval <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))

for (i in 1:length(results)) {
  network = as.numeric(substr(results[[i]]$network, nchar(results[[i]]$network), nchar(results[[i]]$network)))
  d <- data[[network]][, -c(1,2)]
  true = EdgeToAdj(gold[[network]], colnames(d), attr="edge")
  sym.true = isSymmetric(true)
  results[[i]]$inf.net$Weight <- as.numeric(results[[i]]$inf.net$Weight)

  method <- paste("BNFinder", results[[i]]$params, sep="")

  e = NULL
  if (emethod=="netb") {
    if (!sym.true) {
      e <- eval_netb(EdgeToAdj(results[[i]]$inf.net, colnames(d), attr="Weight"),
                     true, sym=FALSE, no.top=no.top) 
    }
    
    e.sym <- eval_netb(EdgeToAdj(results[[i]]$inf.net, colnames(d), attr="Weight"),
                       true, sym=TRUE, no.top=no.top) 
  }
  if (emethod=="mnet") {
    if (!sym.true) {
      e <- eval_minet(EdgeToAdj(results[[i]]$inf.net, colnames(d), attr="Weight"),
                      true, sym=FALSE)
    }
    e.sym <- eval_minet(EdgeToAdj(results[[i]]$inf.net, colnames(d), attr="Weight"),
                        true, sym=TRUE)
  }
  if (!is.null(e)) {
    eval <- rbind(eval, c(method, e[1], e[2], network, results[[i]]$time.sec))
  }
  eval.sym <- rbind(eval.sym, c(method, e.sym[1], e.sym[2], network, results[[i]]$time.sec))
}

methods <- c("FastBMA", "aracne","c3net","clr",
             "Genie3","mrnet",
             "mutrank","mrnetb","pcit")
functions <- c("FastBMA", "aracne.wrap","c3net.wrap","clr.wrap",
               "Genie3.wrap","mrnet.wrap",
               "mutrank.wrap","mrnetb.wrap","pcit.wrap")

for (i in 1:length(data)) {
  network=i
  d <- data[[network]][, -c(1,2)]
  nTimePoints <- length(unique(data[[network]]$time))
  true = EdgeToAdj(gold[[network]], colnames(d), attr="edge")
  sym.true = isSymmetric(true)

  for (j in 1:length(methods)) {  
    method = methods[j]
    print(method)
    start_time <- Sys.time()
    if (method=="FastBMA") {
      inf <- do.call(functions[j], list(d, nTimePoints, priors = reg.prob, known = NULL))
    } else {
      inf <- do.call(functions[j], list(d))
    }
    end_time <- Sys.time()
    time <- as.numeric(end_time - start_time)
    
    sym.inf = isSymmetric(inf)
    print(sym.inf)
    e = NULL  
    if ((sym.true) || (sym.inf)) {
      if (emethod=="netb") {
        e.sym <- eval_netb(inf, true, sym=TRUE, no.top=no.top) 
      }
      if (emethod=="mnet") {
        e.sym <- eval_minet(inf, true, sym=TRUE)
      }
    } else if ((!sym.true)&&(!sym.inf)) {
      if (emethod=="netb") {
        e.sym <- eval_netb(inf, true, sym=TRUE, no.top=no.top) 
        e <- eval_netb(inf, true, sym=FALSE, no.top=no.top) 
      }
      if (emethod=="mnet") {
        e.sym <- eval_minet(inf, true, sym=TRUE)
        e <- eval_minet(inf, true, sym=FALSE)
      }
    }
    
    eval.sym <- rbind(eval.sym, c(method, e.sym[1], e.sym[2], network, time))
    if (!is.null(e)) {
      eval <- rbind(eval, c(method, e[1], e[2], network, time))
    }
  }
}

eval <- eval[-1, ]
eval.sym <- eval.sym[-1, ]
eval <- eval[order(eval$no.net),]
eval.sym <- eval.sym[order(eval.sym$no.net),]

write.table(eval, file = paste("../eval/", data.name, "_eval_",emethod,".tsv", sep=""),
            quote=FALSE, sep="\t", col.names=TRUE, row.names = FALSE)

write.table(eval.sym, file = paste("../eval/", data.name, "_eval_sym_",emethod,".tsv", sep=""),
            quote=FALSE, sep="\t", col.names=TRUE, row.names = FALSE)
