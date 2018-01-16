exprs <- matrix(rexp(9, rate=.1), ncol=3)
colnames(exprs) <- c("G1", "G2", "G3")

exprs <- round(exprs, 3)
cor <- round(cor(exprs), 3)

library(igraph)
g <- graph.adjacency(cor, weighted=TRUE, mode="directed")
plot.igraph(g,vertex.label=V(g)$name,layout=layout.circle,
            vertex.size=37, vertex.label.cex = 1.3, edge.arrow.size=0.5,
            edge.color="black",edge.width=abs(E(g)$weight^2))
