library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
library(minet)
library(GENIE3)
data(dream4)
data(vignette)

Yeast.TFs <- read.table("../data/RegulationTwoColumnTable_Documented_2013927.tsv",
                        check.names = FALSE, header=FALSE, sep=";", quote="")
Yeast.TFs <- levels(unique(Yeast.TFs[,1]))
gene_to_ORF <- select(org.Sc.sgd.db, keys=Yeast.TFs, columns = c("ORF"),
                      keytype="COMMON")

EdgeToAdj <- function(edgeList, colnames, attr) {
  g <- graph.data.frame(edgeList)
  adj <- get.adjacency(g, attr=attr,sparse=FALSE)
  if (length(colnames) > length(colnames(adj))) {
    adj <- expandInfMatrix(adj, colnames)
  }
  
  return(adj)
}

expandInfMatrix <- function(adj, colnames) {
  missing <- colnames[(!(colnames %in% colnames(adj)))]
  #print(length(missing))
  if (length(missing)==0) {
    missing <- colnames[(!(colnames %in% rownames(adj)))]
  }
  adj <- as.data.frame(adj)
  adj[, missing] <- 0
  adj[missing, ] <- 0
  adj <- as.matrix(adj)
  adj <- adj[colnames,colnames]
  return(adj)
}

FastBMA <- function(data, nTimePoints, priors = NULL, known = NULL){
  edges <- networkBMA(data = data, nTimePoints = nTimePoints, prior.prob = priors, known = known,
                      control=fastBMAcontrol(fastgCtrl=fastgControl(optimize=4)))
  g <- graph.data.frame(edges)
  adj <- get.adjacency(g, attr='PostProb',sparse=FALSE)
  
  if (length(colnames(data)) > length(colnames(adj))) {
    adj <- expandInfMatrix(adj, colnames(data))
  }
  
  adj <- adj/max(adj)
  return(adj)
}

eval_minet <- function(inf, gold, sym) {
  if (sym==FALSE) {
    e <- validate(inf, gold)
  } else {
    gold <- pmax(gold,t(gold))
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
emethod <- eval_methods[1]

datasets <- c("brem", "YeastTS", "gnw2000", "gnw2000short")
dataset <- datasets[4]

no.top = 0.2
if (dataset=="brem") {
  data <- t(brem.data)
  gold <- referencePairs
  gold <- gold[which(gold$Regulator %in% colnames(data)),]
  gold <- gold[which(gold$TargetGene %in% colnames(data)),]

  gold$edge <- 1
  gold <- EdgeToAdj(gold, colnames(data), attr="edge")
  sym.true = isSymmetric(gold)
  data.name = "brem"
  methods <- c("aracne","c3net","clr",
               "Genie3", "Genie3.noregs","mrnet",
               "mutrank","mrnetb","pcit")
  functions <- c("aracne.wrap","c3net.wrap","clr.wrap",
                 "GENIE3", "Genie3.wrap","mrnet.wrap",
                 "mutrank.wrap","mrnetb.wrap","pcit.wrap")
}

if (dataset=="YeastTS") {
  data = timeSeries
  nTimePoints <- length(unique(data$time))
  data = data[,-c(1:2)]
  gold = referencePairs
  gold$edge <- 1
  gold = EdgeToAdj(gold, colnames(data), attr="edge")
  sym.true = isSymmetric(gold)
  data.name = dataset
  
  methods <- c("FastBMA","aracne","c3net","clr",
               "Genie3.noregs","mrnet",
               "mutrank","mrnetb","pcit")
  functions <- c("FastBMA","aracne.wrap","c3net.wrap","clr.wrap",
                 "Genie3.wrap","mrnet.wrap",
                 "mutrank.wrap","mrnetb.wrap","pcit.wrap")
}

if ((dataset=="gnw2000")||(dataset=="gnw2000short")) {
  if (dataset=="gnw2000short") {
    print("short")
    data = read.table("../bnf_in/gnw2000short.in", header=TRUE, sep = "\t")
    data = t(data)
  } else {
    data = gnw2000.data
  }
  gold = gnw2000.net
  sym.true = isSymmetric(gold)
  data.name = dataset
  
  methods <- c("aracne","c3net","clr",
               "Genie3", "Genie3.noregs","mrnet",
               "mutrank","mrnetb","pcit")
  functions <- c("aracne.wrap","c3net.wrap","clr.wrap",
                 "GENIE3", "Genie3.wrap","mrnet.wrap",
                 "mutrank.wrap","mrnetb.wrap","pcit.wrap")
}

load(paste("../bnf_out/", data.name, "_results.RData", sep=""))

eval.sym <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))
eval <- data.frame(method=c(1), AUROC=c(1), AUPR=c(1), no.net=c(1), time.sec=c(1))


for (i in 1:length(results)) {
  results[[i]]$inf.net$Weight <- as.numeric(results[[i]]$inf.net$Weight)
  
  method <- paste("BNFinder", results[[i]]$params, sep="")
  print(method)
  e = NULL
  if (emethod=="netb") {
    if (!sym.true) {
      e <- eval_netb(EdgeToAdj(results[[i]]$inf.net, colnames(data), attr="Weight"),
                     gold, sym=FALSE, no.top=no.top) 
    }

    e.sym <- eval_netb(EdgeToAdj(results[[i]]$inf.net, colnames(data), attr="Weight"),
                   gold, sym=TRUE, no.top=no.top) 
  }
  if (emethod=="mnet") {
    if (!sym.true) {
      e <- eval_minet(EdgeToAdj(results[[i]]$inf.net, colnames(data), attr="Weight"),
                    gold, sym=FALSE)
    }
    e.sym <- eval_minet(EdgeToAdj(results[[i]]$inf.net, colnames(data), attr="Weight"),
                    gold, sym=TRUE)
  }
  if (!is.null(e)) {
    eval <- rbind(eval, c(method, e[1], e[2], data.name, results[[i]]$time.sec))
  }
  eval.sym <- rbind(eval.sym, c(method, e.sym[1], e.sym[2], data.name, results[[i]]$time.sec))
}

for (j in 1:length(methods)) {
  method = methods[j]
  print(method)
  start_time <- proc.time()
  if (method=="FastBMA") {
    inf <- do.call(functions[j], list(data, nTimePoints, priors = reg.prob, known = NULL))
  } else if (method=="Genie3") {
    regs <- colnames(data)[which(colnames(data) %in% gene_to_ORF$ORF)]
    inf <- do.call(functions[j], list(t(data), regulators = regs))
    if (length(colnames(data)) > min(length(colnames(inf)), length(rownames(inf)))) {
      inf <- expandInfMatrix(inf, colnames(data))
    }
  } else {
    inf <- do.call(functions[j], list(data))
  }
  end_time <- proc.time()
  time <- (end_time - start_time)[3]

  sym.inf = isSymmetric(inf)
  print(sym.inf)
  e = NULL  
  if ((sym.true) || (sym.inf)) {
    if (emethod=="netb") {
      e.sym <- eval_netb(inf, gold, sym=TRUE, no.top=no.top) 
    }
    if (emethod=="mnet") {
      e.sym <- eval_minet(inf, gold, sym=TRUE)
    }
  } else if ((!sym.true)&&(!sym.inf)) {
      if (emethod=="netb") {
        e.sym <- eval_netb(inf, gold, sym=TRUE, no.top=no.top) 
        e <- eval_netb(inf, gold, sym=FALSE, no.top=no.top) 
      }
      if (emethod=="mnet") {
        e.sym <- eval_minet(inf, gold, sym=TRUE)
        e <- eval_minet(inf, gold, sym=FALSE)
      }
  }

  eval.sym <- rbind(eval.sym, c(method, e.sym[1], e.sym[2], data.name, time))
  if (!is.null(e)) {
    eval <- rbind(eval, c(method, e[1], e[2], data.name, time))
  }
}

eval <- eval[-1, ]
eval.sym <- eval.sym[-1, ]

write.table(eval, file = paste("../eval/", data.name, "_eval_",emethod,".tsv", sep=""),
            quote=FALSE, sep="\t", col.names=TRUE, row.names = FALSE)

write.table(eval.sym, file = paste("../eval/", data.name, "_eval_sym_",emethod,".tsv", sep=""),
            quote=FALSE, sep="\t", col.names=TRUE, row.names = FALSE)
