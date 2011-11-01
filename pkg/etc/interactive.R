ids <- "2009-06-30-15-56-18"
result.dirs <- paste("results/",ids,"/",sep='')
names(result.dirs) <- ids
source("db.R") # for read.res.tables
results <- read.results(result.dirs)
pfiles <- paste(result.dirs,"parameters.txt",sep='/')
params <- combine(lapply(pfiles,read.params))
rownames(params) <- ids
r <- results[[1]]
d <- data.frame(estimate=r$pi.est$MOY,actual=r$pi.sim[,1],r$s[,1:2],
                ppp=r$ppp[r$ppp$POP==1,"PPPVALtot"])
dim(d)
head(d)

library(iplots)
attach(d)
iplot(actual,estimate)
imosaic(factor(s),type)
ihist(ppp)

library(rggobi)
g <- ggobi(d)


## what I would like to write to implement an animation with custom
## interaction that updates the selected locus in each panel
## (currently not possible)
vpl <- function(x,y)viewport(layout.pos.row=x,layout.pos.col=y)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
## let freqs be the data.frame of simulated frequencies for all
## populations, and anc.est be the estimated ancestral allele
## frequencies for each generation
Selector <- makeSelector(locus=1) ## starting selection
Selected.Freqs <- Selector$subset(freqs)
Selected.Est <- Selector$subset(anc.est)
freq.time <- ggplot(Selected.Locus,aes(generation,simulated))+
  geom_line(aes(group=population,colour=colour))+
  geom_vline(aes(time=generation))
print(p,vp=vpl(1,1))
aep <- ggplot(anc.est,aes(simulated,estimated,time=generation))+
  geom_point(aes(colour=type,symbol=type),
             onClick=Selector$select)+
  geom_point(data=Selected.Est,pch=20,cex=2)
pushViewport(vpl(1,2));print(aep,newpage=FALSE);popViewport()
fe <- ggplot(freqs,aes(locus,simulated))+
  geom_point(aes(colour=type,time=generation),
             onClick=Selector$select)+
  geom_vline(data=Selected.Freqs)+
  facet_grid(~type)
## missing components: 1. The magic Selector which makes dynamic
## subsets of data frames. cf mutaframes and the plumbr package?
## 2. the time aesthetic for each value will be a frame of the
## animation. cf QtAnimate and associated web pages............
## 3. the onClick method of geoms. The idea is that you click on the
## plot and then it uses nearest neighbors in pixel space to find the
## nearest plotted data point. Then it calls the specified function
## with that data line, like this fun(locus=1,x=100,y=20).so then we
## filter the displayed results 4. multipanel facets/plots in qt?
## 5. the entire ggplot2 infrastructure on top of qt?
Selector$select <- function(self,...){
  L <- list(...)
  value <- L[[self$select.col]]
  if(is.null(value))value <- L[[1]] ## for fun(1) case
  select.data <- self$data[,self$select.col]
  to.select <- select.data == value
  self$selected <- self$data[to.select,]
}
