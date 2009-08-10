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

