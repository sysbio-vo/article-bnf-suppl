library(ggplot2)
library(grid)
library(dplyr)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

norm = 3600
units = "hours"
dnames = c("sachs_res", "input8_res")
dname = dnames[2]
data <- read.table(paste("../performance/", dname, ".tsv", sep=""))


colnames(data) <- c("algorithm", "cores", 1, 2, 3, 4, 5)
time <- data[3:7]
data <- data[,-c(3:7)]

performance <- data
speedup <- data
efficiency <- data

performance$mean <- apply(time/norm, 1, mean)
performance$se <- apply(time/norm, 1, sd)

sp <- apply(time, 2, function(x) {x[1]/x} )
speedup$mean <- apply(sp, 1, mean)
speedup$se <- apply(sp, 1, sd)

ef <- apply(time, 2, function(x) {x[1]/x/data$cores} )
efficiency$mean <- apply(ef, 1, mean)
efficiency$se <- apply(ef, 1, sd)

pd <- position_dodge(0)
pr <- ggplot(performance %>% arrange(algorithm), aes(x=cores, y=mean, colour=algorithm)) + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=pd, colour="#373738") +
      geom_line(position=pd) + geom_point(aes(shape=algorithm, colour=algorithm)) +
      ylab(paste("Running time,", units)) + ggtitle("Performance comparison") +
      theme(legend.position="bottom", legend.title=element_blank()) +
      scale_x_continuous(breaks=c(data$cores), labels=c(data$cores), minor_breaks = NULL)  

sp <- ggplot(speedup %>% arrange(algorithm), aes(x=cores, y=mean, colour=algorithm)) + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=pd, colour="#373738") +
      geom_line(position=pd) + geom_point(aes(shape=algorithm, colour=algorithm)) +
      ylab("Speedup") + ggtitle("Speedup comparison") +
      theme(legend.position="bottom", legend.title=element_blank()) +
      scale_x_continuous(breaks=c(speedup$cores), labels=c(speedup$cores), minor_breaks = NULL)  

ef <- ggplot(efficiency %>% arrange(algorithm), aes(x=cores, y=mean, colour=algorithm)) + 
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=pd, colour="#373738") +
      geom_line(position=pd) + geom_point(aes(shape=algorithm, colour=algorithm)) +
      ylab("Efficiency") + ggtitle("Efficiency comparison") +
      theme(legend.position="bottom", legend.title=element_blank()) +
      scale_x_continuous(breaks=c(efficiency$cores), labels=c(efficiency$cores), minor_breaks = NULL)  

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(sp)
lwidth <- sum(mylegend$width)

pdf(paste("../plots/", dname, ".pdf", sep=""), height=4, width=12.5)
p4 <- grid.arrange(arrangeGrob(pr + theme(legend.position="none"),
                               sp + theme(legend.position="none"),
                               ef + theme(legend.position="none"),
                               nrow=1),
                   mylegend, heights=c(10, 1), nrow=2)
dev.off()



norm = 60
units = "mins"
dnames = c("r7g7o200_res", "r8g8o10_res", "r8g8o100_res", "r7g7o100_res")
data1 <- read.table(paste("../performance/", dnames[1], ".tsv", sep=""))
data1$dataset <- "7 regulators, 7 targets, 25600 observations"
data2 <- read.table(paste("../performance/", dnames[2], ".tsv", sep=""))
data2$dataset <- "8 regulators, 8 targets, 2560 observations"
data3 <- read.table(paste("../performance/", dnames[3], ".tsv", sep=""))
data3$dataset <- "8 regulators, 8 targets, 25600 observations"
data4 <- read.table(paste("../performance/", dnames[4], ".tsv", sep=""))
data4$dataset <- "7 regulators, 7 targets, 12800 observations"


data <- rbind(data1, data2, data3, data4)
colnames(data) <- c("algorithm", "cores", 1, 2, 3, 4, 5, "dataset")
data <- data[which(data$algorithm=="hybrid"),]

time <- data[3:7]
data <- data[,-c(3:7)]
performance <- data
performance$mean <- apply(time/norm, 1, mean)
performance$se <- apply(time/norm, 1, sd)
performance$cores <- as.factor(performance$cores)
performance$dataset <- factor(performance$dataset,levels(factor(performance$dataset))[c(3,2,1,4)])
performance <- performance[which(performance$cores %in% c(1, 2)),]

pd <- position_dodge(width = 0.9)

pl <- ggplot(performance, aes(x=cores, y=mean, fill=dataset)) +
  geom_bar(stat="identity", position=pd) +
  geom_text(aes(label=round(mean, 2)), vjust=-0.6, color="#0e0e0f",
            position = position_dodge(0.9), size=3.5) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="#373738",
                position=pd, width = 0.25, alpha=0.6)+
  #scale_fill_brewer(palette="Set2") +
  scale_fill_manual(values=c("#33D47F","#FDF06E","#4A9AF6","#FB6099")) +
  ylab("Running time, mins")

pl
library(cowplot)
theme_set(theme_grey())
save_plot("../plots/datasets_perf_comparison.pdf", pl, base_height=5, base_width=12, nrow=1)
