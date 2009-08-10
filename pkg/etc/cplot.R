cc <- results[[1]]$c
clong <- reshape(cc,dir="long",sep="",varying=colnames(cc),timevar="pop")
names(clong)[2] <- "c"
s <- results[[1]]$s
s$id <- rownames(s)
ltab <- merge(s,clong)
densityplot(~c,ltab,groups=type,main="Estimate of differentiation does not depend on selection type")
