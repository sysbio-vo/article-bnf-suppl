## Testing netbenchmark

library(devtools)
dev_mode(on=T)
install_github("paubellot/netbenchmark", ref="master")
library(netbenchmark)

ngenes <- dim(syntren300.data)[2]
nlinks <- ngenes^2-ngenes
no.topedges=100
no.edges <- round(nlinks*no.topedges/100)

# True net versus true net
r1 <- evaluate(syntren300.net,syntren300.net,extend=no.edges,sym=FALSE)
aupr(r1)
auroc(r1)
r2 <- evaluate(syntren300.net,syntren300.net,extend=0,sym=FALSE)
aupr(r2)
auroc(r2)

# True net with 0.9 instead of zeros versus true net
syntren300.net.nozeros <- syntren300.net
syntren300.net.nozeros[syntren300.net.nozeros==0] <- 0.9
syntren300.net.nozeros[4, ] <- 0

r3 <- evaluate(syntren300.net.nozeros,syntren300.net,extend=0,sym=FALSE)
aupr(r3)
auroc(r3)
r4 <- evaluate(syntren300.net.nozeros,syntren300.net,extend=no.edges,sym=FALSE)
aupr(r4)
auroc(r4)

dev_mode(on=F)

## Testing minet

library(minet)

# True net versus true net
e1 <- validate(syntren300.net,syntren300.net)
auc.pr(e1)
auc.roc(e1)

# True net with 0.9 instead of zeros versus true net
e2 <- validate(syntren300.net.nozeros,syntren300.net)
auc.pr(e2)
auc.roc(e2)

## Testing networkBMA

library(networkBMA)
library(igraph)

# True net versus true net
g <- graph.adjacency(syntren300.net, weighted=TRUE, mode="undirected")
true <- as.data.frame(get.edgelist(g))
true$weight <- E(g)$weight

contingencyTables <- contabs(network = true, reference = true[,1:2],
                             size = no.edges)
prc(contingencyTables, plotit = FALSE)
roc(contingencyTables, plotit = FALSE)

# True net with 0.9 instead of zeros versus true net
g <- graph.adjacency(syntren300.net.nozeros, weighted=TRUE, mode="undirected")
inf <- as.data.frame(get.edgelist(g))
inf$weight <- E(g)$weight

contingencyTables <- contabs(network = inf, reference = true[,1:2],
                             size = no.edges)
prc(contingencyTables, plotit = FALSE)
roc(contingencyTables, plotit = FALSE)
