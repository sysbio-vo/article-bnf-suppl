library(netbenchmark)
library(networkBMA)
library(igraph)
library(org.Sc.sgd.db)
library(abind)
data(dream4)
data(vignette)


BNFinder <- function(data, data.name, regulators=NULL, priors=NULL, lim=0, sub=0, k=0, path=""){
  in.name <- paste(getwd(), "/../bnf_in/", data.name, ".in", sep="")
  write.table(t(data), in.name, sep="\t", quote=FALSE)

  if (!is.null(regulators)) {
    fConn <- file(in.name, 'r+')
    Lines <- readLines(fConn) 
    firstLine <- paste("#regulators", paste(regulators, collapse = " "))
    writeLines(c(firstLine, Lines), con = fConn) 
    close(fConn)
  }

  if (!is.null(priors)) {
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

  if (k>0) {
    time <- read.table(paste(getwd(), "/time", k, ".txt", sep=""))
  }
  return(list(time, res))
}


## Yeast
# Read Yeast TFs
Yeast.TFs <- read.table("../data/RegulationTwoColumnTable_Documented_2013927.tsv",
                        check.names = FALSE, header=FALSE, sep=";", quote="")
Yeast.TFs <- levels(unique(Yeast.TFs[,1]))
gene_to_ORF <- select(org.Sc.sgd.db, keys=Yeast.TFs, columns = c("ORF"),
                      keytype="COMMON")

regs <- colnames(gnw2000.data)[which(colnames(gnw2000.data) %in% gene_to_ORF$ORF)]


## Dream4 TS

params <- data.frame(lim=c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3),
                     sub = c(0, 5, 10, 20, 30, 0, 5, 10, 20, 30, 0, 5, 10, 20, 30, 0, 5, 10, 20, 30))

params <- data.frame(lim=c(1), sub=c(5) )

data = dream4ts10
data.name = "dream4ts10"
path = paste(getwd(), "/../bnfinder/", sep="")

results <- c()
for (network in 1:length(data)) {
  print(paste("Processing network", network))
  dat <- data[[network]][, -c(1:2)]
  for (i in 1:nrow(params)) {
    net.name <- paste(data.name, "n", network, sep="")
    result <- BNFinder(dat, net.name, lim=params[i, 1], sub=params[i, 2], k=2, path=path)
    results <- c(results,
                 list(list(network=net.name, time.sec=as.numeric(result[[1]]), inf.net=result[[2]])))
  }
}  

save(results, file="results.RData")