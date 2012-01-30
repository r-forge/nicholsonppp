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

## example of how to do animation and interactivity at the same time.
library(qtbase)
b <- Qt$QPushButton("something")
qconnect(b,"pressed",function()print("foo"))
b$show()
a <- Qt$QPropertyAnimation(b,"geometry")
a$setDuration(20000)
width <- 100
height <- 30
x <- 830
y <- 200
a$setStartValue(Qt$QRect(x,y,width,height))
a$setEndValue(Qt$QRect(x+250,y+250,width+100,height+100))
a$start()
a


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

## it would have been nice to write this to implement the
## breakpoint/copies annotator that took 400 lines of python code
ggplot()+
  aes(position,logratio)+
  geom_point(aes(colour=copies,fill=annotation),data=profiles,
             onClick=cycleAnnotation)+
  geom_tallrect(aes(fill=annotation),ymin=0.5,ymax=1,data=copies,
                onClick=cycleDeleteAnn)+
  geom_tallrect(aes(fill=annotation),ymin=0,ymax=0.5,data=breakpoints,
                onClick=cycleDeleteAnn)+
  facet_grid(profile_id~chromosome)+
  interaction(ymin=0.5,ymax=1,onDrag=TallRectangle,
              onRelease=copies$add_row)+
  interaction(ymin=0,ymax=0.5,onDrag=TallRectangle,
              onRelease=breakpoints$add_row)+
  interaction(rightClick=exchangeProfile)

## again, with the annotated copy number profiles, it would be nice to
## have 1 plot of error curves linked to a plot of the smoothing
## model... i.e. click the profiles to add annotations, which get used
## to construct the error curves. Click the error curve to update the
## smoothing model, superimposed on the profiles.
## dimensions: profiles, chromosomes, parameters, positions
DataSet <- proto(.,{
  select
  insert
  update
  delete
})
Geom <- proto(.,{
  on_click
})
## parameters is a table with columns id, parameter, error
selected_parameter <- Selector(parameters,id)
makeAnnotation <- list(dragging=function(x.begin,x.now,...){
  resize_tallrect(x.begin,x.now)
},mouseUp=function(x.begin,x.now,chromosome,profile_id,...){
  beg <- to_position(x.begin)
  end <- to_position(x.now)
  m <- min(beg,end)
  M <- max(beg,end)
  annotations$insert(profile_id,chromosome,m,M,"breakpoint")
  ## insert function should look for attached plots and send them
  ## update signals
  region <- smooth$select(profile_id,chromosome,m,M)
  for(p in unique(region$parameter)){
    region.smooth <- subset(region,parameter==p)
    ## store agreement for this region somewhere
    parameters$update() ## update for each parameter
  }
})
cycleAnnotation <- function(geom,...){
  geom$table ## should be a link to the dataset
  current_ann <- geom$get("annotation")
  if(current_ann%in%names(next_annotation)){
    geom$update(next_annotation[current_ann])
  }else{
    geom$delete()
  }
}
profile_plot <- ggplot()+
  facet_grid(profile_id~chromosome)+
  geom_point(aes(position,logratio),data=probes)+
  geom_line(aes(position,smooth),data=smooth)+
  geom_tallrect(aes(min=min,max=max),data=annotations,
                interactions=list(click=cycleAnnotation))+
  interactions(drag=makeAnnotation)
error_plot <- ggplot()+
  geom_line(aes(log10(parameter),error),data=parameters)+
  geom_vline(aes(xintercept=parameter),data=selected_parameter)+
  interactions(click=selected_parameter$select)
## then plot them both... start with no model hilited, no annotations,
## flat error curve. then drag the profile to add annotations, and the
## error curve should get updated. makeAnnotation must send a signal
## to update the error column to the parameters data table, which is
## now linked to the error curve geom_line. Then click the error curve
## to set a vertical line, highlighting a model parameter, and sending
## a signal to the smooth data set, which in turn updates the drawn
## lines...
