library(ggplot2)
library(fgsea)
plotEnrichment2<-function (pathway, stats, gseaParam=1, ticksSize=0.2,linecol=NULL,colticks="#00000099") 
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats=pathway, 
                          returnAllExtremes=TRUE)
  # Set line color
  if(is.null(linecol)){
    if(gseaRes$res>=0){
      linecol<-"red3"
    } else {
      linecol<-"navy"
    }
  }
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x=y=NULL
  g<-ggplot(
    toPlot,
    aes(x=x,y=y)) + 
    geom_point(color=linecol,size=0.1) +
    geom_hline(yintercept=max(tops),colour="red",linetype="dashed") +
    geom_hline(yintercept=min(bottoms),colour="red",linetype="dashed") +
    geom_hline(yintercept=0,colour="black") +
    geom_line(color=linecol) + theme_bw() +
    geom_segment(data=data.frame(x=pathway),mapping=aes(x=x,y=-diff/2,xend=x,yend=diff/2),size=ticksSize,col=colticks) + 
    theme(panel.border=element_blank(), panel.grid.minor=element_blank()) + 
    labs(x="rank", y="enrichment score")
  g
  
  # TODO in the future: a lattice way to plot GSEA (one-way and two-ways)
}
