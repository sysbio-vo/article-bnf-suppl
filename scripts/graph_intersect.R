# BNF output comparison
library(igraph)


gnw1 <- read.table("../bnf_out/gnw2000L1I40.out", header=F)
gnw2 <- read.table("../bnf_out/gnw2000L2I40.out", header=F)
gnw3 <- read.table("../bnf_out/gnw2000L1I20.out", header=F)

gnw1 <- gnw1[order(gnw1$V2, decreasing = T),]
gnw2 <- gnw2[order(gnw2$V2, decreasing = T),]

g1 <- graph_from_data_frame(gnw1[,c(1, 3)])
g2 <- graph_from_data_frame(gnw2[,c(1, 3)])

intersect <- g1 %s% g2
length(E(intersect))
