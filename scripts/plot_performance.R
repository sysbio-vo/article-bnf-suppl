library(ggplot2)
library(grid)
library(dplyr)
library(gridExtra)
library(reshape2)
library(RColorBrewer)

norm = 60
dname = "sachs"
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
      ylab("Performance, mins") + ggtitle("Performance comparison") +
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
