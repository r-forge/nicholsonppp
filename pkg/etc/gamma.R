source("db.R")
scvals <- c(0.1,0.5,1)
shvals <- c(.1,.05,.01,.005,.001)
x <- seq(0,max(sapply(scvals,function(v)v*shvals))+1,l=1e3)
d <- NULL
for(scale in scvals)for(shape in shvals)
  d <- rbind(d,
             data.frame(value=x,
                        density=dgamma(x,shape,scale=scale),
                        scale,
                        shape))
panel.densitycompare <- function(x,y,subscripts,...){
  p <- d[subscripts[1],c("shape","scale")]
  grid.lines(p$shape*p$scale,
             c(0,1),
             default.units="native",
             gp=gpar(col='grey'))
  grid.lines(qgamma(0.99,shape=p$shape,scale=p$scale),
             c(0,1),
             default.units="native",
             gp=gpar(col='red'))
  cv <- 1
  prob <- 1-pgamma(cv,shape=p$shape,scale=p$scale)
  pos <- 0.98
  figs <- 10
  grid.text(format(round(prob,digits=figs),ns=figs,scientific=F),pos,pos,
            just=c("right","top"),
            gp=gpar(col='lightblue')
            )
  panel.xyplot(x,y,subscripts=subscripts,...)
  y <- c(0,y[x>cv],0)
  x <- x[x>cv]
  x <- c(x[1],x,x[length(x)])
  grid.polygon(x,y,default.units="native",gp=gpar(fill='lightblue',lty=0))
}
mypng('gamma-densities',
xyplot(density~value|scale+shape,d,
            type='l',
            strip=strip.custom(strip.names=c(T,T),strip.levels=c(T,T)),
            ylim=c(0,0.2),
            xlim=c(0,max(x)),
            panel=panel.densitycompare,
            main='Gamma densities',
            sub="99%tile in red, mean in grey, P[X>1] shaded and displayed in blue"),
      paper='a4'
      )

