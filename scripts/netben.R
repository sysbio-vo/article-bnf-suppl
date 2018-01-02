library(netbenchmark)


# Define a wrapper function
Spearmancor <- function(data, regulators){
  regulators = 1
  cor(data,method="spearman")
}
## Not run:
# Register it to all.fast methods
RegisterWrapper("Spearmancor")
# Register it to all methods
RegisterWrapper("Spearmancor", all.fast=FALSE)



top20 <- netbenchmark(seed=42, datasets.num=1, experiments = 2001)
top20.fast.list[[1]]

View(rogers1000.data)
View(rogers1000.net)

View(gnw2000.net)

Yeast.TFs <- read.table("RegulationTwoColumnTable_Documented_2013927.tsv",
                          check.names = FALSE, header=FALSE, sep=";", quote="")

Yeast.TFs <- levels(unique(Yeast.TFs[,1]))

library(igraph)
g <- graph.adjacency(gnw2000.net)
edg <- get.edgelist(g)


library(org.Sc.sgd.db)
gene_to_ORF <- select(org.Sc.sgd.db, keys=Yeast.TFs, columns = c("ORF"),
                      keytype="COMMON")

genes <- unique(edg[,1])

length(which(gene_to_ORF$ORF %in% genes))

library(gsubfn)
pat <- read.pattern("TFConsensusList_20130918.transfac", pattern="^NA + (\\S+)")


Lines <- "this is the first field 1 2
more text 3 4
"
pat <- "^(.*) +(\\S+) +(\\S+)$"
read.pattern(text = Lines, pattern = pat, as.is = TRUE)



library(networkBMA)
data(dream4)
dream4ts10[[1]][,-(1:2)]
dream4ts10[[1]][,3]

data(vignette)
ScanBMA()
networkBMA()

head(brem.data)

network=1
nTimePoints <- length(unique(dream4ts10[[network]]$time))
edges1ts10 <- networkBMA( data = dream4ts10[[network]][,-(1:2)], 
                          nTimePoints = nTimePoints, prior.prob = 0.01)
summary(edges1ts10)

ctables <- contabs.netwBMA(edges1ts10, dream4gold10[[network]],
                           reg.known, thresh=c(.5,.75,.9))