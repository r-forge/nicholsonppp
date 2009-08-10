##ids <- "2009-06-23-14-57-07"
result.dirs <- paste("results/",ids,"/",sep='')
names(result.dirs) <- ids
source("db.R") # for read.res.tables
results <- read.results(result.dirs)
pfiles <- paste(result.dirs,"parameters.txt",sep='/')
params <- combine(lapply(pfiles,read.params))
rownames(params) <- ids

r <- results[[1]]
a <- cbind(est.naive=rowMeans(r$alpha.sim),
           est.Nicholson=r$pi.est[,2],
           r$s,
           actual=r$pi.sim[,1])
a$id <- 1:nrow(a)
head(a)
b <- reshape(a,direction="long",varying=1:2,timevar="method",idvar="locus")
source("db.R")
subdir <- paste("results/",ids,"/",sep="")
mypng("nicholson-naive",
      pi.xyplot(est~actual|method,cbind(b,fixnames(params)),alpha=1,
                main="Naive estimation does not agree with Nicholson model estimates"),
      subdir)

pos <- a[a$type=="POS",]
d <- abs(pos$actual-0.5)
r$alpha.sim[pos[d==min(d),'id'],]

## density estimation for cutoff values
ppp <- r$ppp
dens <- density(ppp$PPPVALtot)
minima <- diff(sign(diff(dens$y)))>0
minx <- dens$x[minima]
pos.cutoff <- minx[1]
bal.cutoff <- minx[length(minx)]
loci <- cbind(r$s,ppp=ppp$PPPVALtot[ppp$POP==1])[,-c(1,3)]
loci$pred <- "NEU"
loci$pred[loci$ppp>bal.cutoff] <- "BAL"
loci$pred[loci$ppp<pos.cutoff] <- "POS"
loci$error <- loci$pred!=loci$type
loci$summary <- paste(loci$pred,loci$type,sep="|")
summary(factor(loci$summary[loci$error]))
table(loci[,c("type","pred")]) ## contingency table
error.rate <- sum(loci$error)/nrow(loci) ## error rate

pdfname <- paste(subdir,"ppp-classify.pdf",sep="")
pdf(pdfname,paper="a4r",w=0,h=0)
plot(dens,main="Classifying locus type based on PPP-values")
mtext(paste("error rate: ",round(error.rate*100,digits=2),"%",sep=""))
abline(v=c(pos.cutoff,bal.cutoff))
ymax <- max(dens$y)
text(pos.cutoff/2,ymax,"POS")
text((bal.cutoff+1)/2,ymax,"BAL")
text((pos.cutoff+bal.cutoff)/2,ymax,"NEU")
dev.off()
system(paste("xpdf",pdfname))


