library(ggplot2)
library(grid)
library(dplyr)
library(gridExtra)
library(reshape2)
library(RColorBrewer)


mname = "mnet"
sym = ""

bnf = "BNFinderL2I30"
dname = "dream4ts10"
data <- read.table(paste("../eval/", dname,"_eval_", sym, mname, ".tsv", sep=""), header=T)

#b <- data[which(data$method==bnf),]
#data <- data[-grep("BNFinder", data$method),]
data <- data[grep("BNFinder", data$method),]
#data <- rbind(data, b)
data <- data[order(data$no.net),]

data <- data %>% group_by(no.net) %>%
        mutate(roc.rank = order(order(AUROC, decreasing=T)))

data <- data %>% group_by(no.net) %>%
        mutate(pr.rank = order(order(AUPR, decreasing=T)))

cols <- colorRampPalette(brewer.pal(9, "Set1"))
myPal <- cols(length(unique(data$method)))

pl.roc <- ggplot(data, aes(x = no.net, y = roc.rank, colour = method)) +
          geom_line(size = 2, alpha = 0.8) +
          geom_point(size = 2.3) +
          geom_point(color = "#FFFFFF", alpha = .8, size = .4) +
          geom_text(data = subset(data, no.net==1), aes(label = method, x=0.8), hjust = 1, size = 4, colour="#888888") +
          scale_x_continuous(breaks = c(1:length(unique(data$no.net))), limits = c(0,length(unique(data$no.net)))) +
          theme(legend.position = "none") +
          labs(x = "Network", y = "Rank") +
          theme(panel.grid.major.y = element_blank()) +
          theme(panel.grid.minor = element_blank()) +
          theme(panel.grid.major.x = element_line(color = "#F3F3F3")) +
          scale_colour_manual(values = myPal) +
          scale_y_reverse(labels=c(1:length(unique(data$method))), breaks=c(1:length(unique(data$method)))) +
          ggtitle(expression(paste(bold("A."), "   AUROC ranking")))
pl.roc
pl.pr <- ggplot(data, aes(x = no.net, y = pr.rank, colour = method)) +
          geom_line(size = 2, alpha = 0.8) +
          geom_point(size = 2.3) +
          geom_point(color = "#FFFFFF", alpha = .8, size = .4) +
          geom_text(data = subset(data, no.net==1), aes(label = method, x=0.8), hjust = 1, size = 4, colour="#888888") +
          scale_x_continuous(breaks = c(1:length(unique(data$no.net))), limits = c(0,length(unique(data$no.net)))) +
          theme(legend.position = "none") +
          labs(x = "Network", y = "Rank") +
          theme(panel.grid.major.y = element_blank()) +
          theme(panel.grid.minor = element_blank()) +
          theme(panel.grid.major.x = element_line(color = "#F3F3F3")) +
          scale_colour_manual(values = myPal) +
          scale_y_reverse(labels=c(1:length(unique(data$method))), breaks=c(1:length(unique(data$method)))) +
          ggtitle(expression(paste(bold("B."), "   AUPR ranking")))

pl.pr
pdf(paste("../plots/", dname,"_", mname,"_", sym, "accuracy.pdf", sep=""), height=4, width=12.5)
p <- grid.arrange(arrangeGrob(pl.roc, pl.pr, nrow=1))
dev.off()


mname = "netb"
sym = "sym_"
bnf = "BNFinderL2I30"
dname = "gnw"
gnw2000 <- read.table(paste("../eval/gnw2000_eval_", sym, mname, ".tsv", sep=""), header=T)
gnw2000short <- read.table(paste("../eval/gnw2000short_eval_", sym, mname, ".tsv", sep=""), header=T)

data <- rbind(gnw2000short, gnw2000)

b <- data[which(data$method==bnf),]
data <- data[-grep("BNFinder", data$method),]
data <- rbind(data, b)
data <- data[order(data$no.net),]

pl.roc <- ggplot(data, aes(x=method, y=AUROC, fill=no.net)) +
          geom_bar(stat="identity", position=position_dodge()) +
          geom_text(aes(label=round(AUROC, 2)), vjust=1.6, color="#0e0e0f",
                position = position_dodge(0.9), size=3.5)+
          scale_fill_manual(values=c("#33D47F","#FB6099")) +
          ggtitle(expression(paste(bold("A."), "   Evaluation by AUROC")))

pl.pr <- ggplot(data, aes(x=method, y=AUPR, fill=no.net)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=round(AUPR, 2)), vjust=1.6, color="#0e0e0f",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values=c("#33D47F","#FB6099")) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  ggtitle(expression(paste(bold("A."), "   Evaluation by AUPR")))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pl.pr)

pdf(paste("../plots/", dname,"_", mname,"_", sym, "accuracy.pdf", sep=""), height=5, width=16)
p <- grid.arrange(arrangeGrob(pl.roc + theme(legend.position="none"),
                              pl.pr + theme(legend.position="none"), nrow=1),
                  mylegend, heights=c(10, 1), nrow=2)
dev.off()



mname = "mnet"
bnf = "BNFinderL2I30"
Yeast = F
if (Yeast) {
  sym = "sym_"
  dname = "YeastTS"
  YeastTS.sym <- read.table(paste("../eval/YeastTS_eval_", sym, mname, ".tsv", sep=""), header=T)
  sym = ""
  YeastTS <- read.table(paste("../eval/YeastTS_eval_", sym, mname, ".tsv", sep=""), header=T)
  YeastTS.sym$eval = "symmetrical"
  YeastTS$eval = "non-symmetrical"
  data <- rbind(YeastTS, YeastTS.sym)
} else {
  dname = "brem"
  sym = "sym_"
  brem.sym <- read.table(paste("../eval/brem_eval_", sym, mname, ".tsv", sep=""), header=T)
  sym = ""
  brem <- read.table(paste("../eval/brem_eval_", sym, mname, ".tsv", sep=""), header=T)
  brem.sym$eval = "symmetrical"
  brem$eval = "non-symmetrical"
  data <- rbind(brem, brem.sym)
}
b <- data[which(data$method==bnf),]
data <- data[-grep("BNFinder", data$method),]
data <- rbind(data, b)
data <- data[order(data$no.net),]

pl.roc <- ggplot(data, aes(x=method, y=AUROC, fill=eval)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=round(AUROC, 2)), vjust=1.6, color="#0e0e0f",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values=c("#FDF06E","#4A9AF6")) +
  ggtitle(expression(paste(bold("A."), "   Evaluation by AUROC")))

pl.pr <- ggplot(data, aes(x=method, y=AUPR, fill=eval)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=round(AUPR, 2)), vjust=1.6, color="#0e0e0f",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(values=c("#FDF06E","#4A9AF6")) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  ggtitle(expression(paste(bold("B."), "   Evaluation by AUPR")))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pl.pr)

pdf(paste("../plots/", dname,"_", mname,"_accuracy.pdf", sep=""), height=5, width=16)
p <- grid.arrange(arrangeGrob(pl.roc + theme(legend.position="none"),
                              pl.pr + theme(legend.position="none"), nrow=1),
                  mylegend, heights=c(10, 1), nrow=2)
dev.off()



mname = "netb"
sym = "sym_"
gnw2000 <- read.table(paste("../eval/gnw2000_eval_", sym, mname, ".tsv", sep=""), header=T)
gnw2000short <- read.table(paste("../eval/gnw2000short_eval_", sym, mname, ".tsv", sep=""), header=T)
dream100.sym <- read.table(paste("../eval/dream4ts100_eval_", sym, mname, ".tsv", sep=""), header=T)
dream10.sym <- read.table(paste("../eval/dream4ts10_eval_", sym, mname, ".tsv", sep=""), header=T)
YeastTS.sym <- read.table(paste("../eval/YeastTS_eval_", sym, mname, ".tsv", sep=""), header=T)
brem.sym <- read.table(paste("../eval/brem_eval_", sym, mname, ".tsv", sep=""), header=T)
sym=""
brem <- read.table(paste("../eval/brem_eval_", sym, mname, ".tsv", sep=""), header=T)
dream100 <- read.table(paste("../eval/dream4ts100_eval_", sym, mname, ".tsv", sep=""), header=T)
dream10 <- read.table(paste("../eval/dream4ts10_eval_", sym, mname, ".tsv", sep=""), header=T)
YeastTS <- read.table(paste("../eval/YeastTS_eval_", sym, mname, ".tsv", sep=""), header=T)

#YeastTS.sym$no.net <- "YeastTS.symmetrical"
#brem.sym$no.net <- "brem.symmetrical"

#dname="brem"
#data <- rbind(brem, brem.sym)

#dname="gnw"
#data <- rbind(gnw2000, gnw2000short)

#dname="dream10"
#data <- rbind(dream10)

#dname="YeastTS"
#data <- rbind(YeastTS, YeastTS.sym)

data <- data[grep("BNFinder", data$method),]
data$limit <- as.numeric(substr(data$method, 10, 10))
data$subopt <- as.numeric(substr(data$method, 12, 14))

data <- data[which(data$limit != 0),]
data[which(data$subopt == 0),]$subopt <- 1

data$limit <- as.factor(data$limit)
data$subopt <- as.factor(data$subopt)

pl <- ggplot(data, aes(x=AUROC, y=AUPR, fill=subopt, colour=subopt)) +
    geom_point(size=4, alpha=0.8, aes(shape = limit)) +
    facet_grid(no.net ~ .) +
    scale_shape_manual(values = c(21, 22, 24))
    
mylegend<-g_legend(pl)
pl <- ggplot(data, aes(x=AUROC, y=AUPR, fill=subopt, colour=subopt)) +
  geom_point(size=4, alpha=0.8, colour="black", aes(shape = limit)) +
  facet_grid(no.net ~ .) +
  scale_shape_manual(values = c(21, 22, 24)) + theme(legend.position="none")

pdf(paste("../plots/", dname,"_", mname,"_arguments.pdf", sep=""), height=5, width=7)
p <- grid.arrange(arrangeGrob(pl, nrow=1),
                  mylegend, widths=c(10, 3), ncol=2, nrow=1)
dev.off()
