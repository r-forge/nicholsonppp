CEX.LEG <- 1
##debug(fixation.endpoints)
display.params=c("populations",'loci.per.s.value','popsize','generations','p.neutral','fixed')
deduce.param.label <- function(lt,imp=display.params){
  L <- sapply(imp,function(cn)if(cn%in%names(lt))
              {tmp <- unique(lt[,cn]);if(is.numeric(tmp))tmp <- round(tmp,digits=2);tmp} else NULL)
  lab <- paste(sapply(names(L),function(N)paste(N,
                                         ": ",
                                         paste(L[[N]],collapse=", "),
                                         "  ",
                                         sep='')),
        collapse='')
  if(lab=="")return(NULL)
  else return(lab)
}
##debug(deduce.param.label)
ltable <- function(N,lt,fixed){
  get.actual <- function(LT)unlist(LT$pi.sim)
  actual <- get.actual(lt)
  if(is.null(actual))actual <- get.actual(results[[N]])[as.integer(gsub("[^0-9]+","",lt$s$old))]
  if(is.null(lt$pi.est))return(NULL)
  data.frame(populations=params[N,'Npop'],
             generations=params[N,'Ngeneration'],
             popsize=params[N,'Taille.POP'],
             loci=params[N,'N.s.rep'],
             fixed=factor(fixed),
             p.neutral=params[N,'p.neutral'],
             actual,
             estimated=lt$pi.est$MOY,
             ppp=lt$ppp$PPPVALtot,
             s=lt$s[,1],
             type=lt$s$type,
             row.names=NULL,
             check.rows=F)
}
##debug(ltable)
make.lattice.table <- function(N){
  df <- ltable(N,results[[N]],'all')
  for(sn in names(results[[N]]$subsets))
    df <- rbind(df,ltable(N,results[[N]]$subsets[[sn]],sn))
  print(N)
  df
}
corr.panel <- function(x,y,...){
  grid.text(paste("N=",length(x),", r=",
                  format(cor(x,y),digits=2,nsmall=2),sep=""),
            0.95,0.05,just=c('right','bottom'))
  panel.xyplot(x,y,...)
}
latticeplot <- function(f,d,FUN,alpha=0.2,omit=NULL,groups=NULL,...){
  fres <- sapply(all.vars(f),function(v)length(grep(v,gsub('^[^|]*','',attr(terms(f),'term.labels')))))
  to.display <- names(d)[names(d)%in% display.params&!(names(d)%in%names(fres)[fres==1])]
  if(is.null(omit))if(!is.null(substitute(groups))){
    ch <- deparse(substitute(groups))
    r <- regexpr('[^()]*[)]',ch,perl=T)
    omit <- if(r==-1)ch else substr(ch,r,r-2+attr(r,'match.length'))
  }
  subset <- -which(to.display%in%omit)
  if(!is.null(omit)&&length(subset))to.display <- to.display[subset]
  mf <- match.call()
  fi <- which(names(mf)=="FUN")
  mf[[1]] <- mf[[fi]]
  mf <- mf[-fi]
  names(mf)[2:3] <- ""
  mf$sub <- deduce.param.label(d,to.display)
  mf$alpha <- alpha
  if(!"strip"%in%names(mf))
    mf$strip <- strip.custom(strip.names=c(T,T),strip.levels=c(T,T))
  eval.parent(mf)
}
##debug(latticeplot)
pi.xyplot <- function(f,d,panel=corr.panel,aspect=1,cols=NULL,...){
  latticeplot(f,d,xyplot,
              groups=type,
              par.settings=list(superpose.symbol=list(
                                  pch=c("B","N","P"),
                                  cex=1.2,
                                  col=cols)),
              panel=panel,
              ylim=c(0,1),
              xlim=c(0,1),
              aspect=aspect,
              ...)
}
## actually creates a PNG for presentations and PDFs for printout
mypng <- function(title,P,subdir='',paper='a4r',show='pdf'){
  pre <- paste(format(Sys.time(),"%Y-%m-%d"),title,sep='-')
  pre <- paste(subdir,pre,sep='')
  f <- paste(pre,'.png',sep='')
  png(f,h=1024*2,w=1280*2,pointsize=100)
  trellis.par.set(fontsize=list(text=40,points=25))
  print(P)
  dev.off()
  
  pdfname <- paste(pre,'.pdf',sep='')
  pdf(pdfname,paper=paper,w=0,h=0)
  print(P)
  dev.off()
  cmd <- paste("eog",f,"&")
  ##print(cmd)
  if(show=='png')system(cmd)
  if(show=='pdf')system(paste("xpdf",pdfname,"&"))
  pdfname
}
fixnames <- function(params){
  names(params) <- rename(names(params),
                          c(N.s.rep='loci',
                            Ngeneration='generations',
                            Taille.POP='popsize',
                            Npop='populations'))
  params
}
fixation.plot <- function(N){
  com <- compare.alleles(N)
  params2 <- fixnames(params)
  params2$id <- rownames(params2)
  com <- merge(com,params2,'id')
  nonzero <- !is.nan(com$prop)
  nonall <- !as.character(com$fixed)=='all'
  mypng('fixation-selection',
        latticeplot(prop~S,
                    com[nonzero&nonall,],
                    xyplot,
                    alpha=1,
                    groups=factor(fixed),
                    ylab='proportion of loci "not fixed"',
                    xlab='selection type and coefficient',
                    main='Loci fixation depends on selection type and coefficient',
                    scales=list(
                      x=list(rot=90,cex=0.65,fontfamily='mono'),
                      y=list(axs='i',rot=0,at=c(0,1,0.25,0.5,0.75))),
                    auto.key=list(space='right'),
                    panel=function(...){grid.lines(c(0.5,0.5),c(0,1),gp=gpar(col='grey'));panel.xyplot(...)},
                    ylim=c(0,1)),
        paste("results/",N,"/",sep='')
        )
}

##debug(make.lattice.table)
combine <- function(L){
  for(i in 1:length(L)){
    if(i==1)ft <- L[[1]]
    else ft <- rbind(ft,L[[i]])
  }
  ft
}
##debug(combine)
read.params <- function(f){
  x <- read.table(f)
  if(!"p.neutral"%in%names(x))x$p.neutral <- 1/3
  x
}
read.s <- function(sf){
  s <- read.table(sf)
  if(ncol(s)==1){
    s$type <- factor(gsub('[0-9]','',rownames(s)))
    colnames(s)[1] <- 's'
    s$oldname <- rownames(s)
    rownames(s) <- NULL
  }
  s
}
## read all data into master data index
read.res.tables <-
  function(d,
           fnames=c(alpha.est="summary_alpha.out",
             alpha.sim="RES.F.fin.txt",
             pi.sim="F.init.txt",
             pi.est="summary_pi.out",
             s="S.LOC.txt",
             c="res_c.out",
             ppp="PPPval.out")){
  files <- paste(d,fnames,sep='')
  names(files) <- names(fnames)
  x <- lapply(files,read.nicholson.outfile)
  x
}
##debug(read.res.tables)


## functions for calculation of loci exclusion
allcols <- function(d,v=NULL,count=ncol(d)){
  if(mode(v)=="function")res <- v(d)
  else res <- matrix(unlist(d)%in%v,nrow=nrow(d),ncol=ncol(d))
  rowSums(res)==count
}
onecol <- function(d,v)allcols(d,v,0)
maf.limit <- 0.01
lo <- function(d)d<maf.limit
hi <- function(d)d>(1-maf.limit)
## order is important here, it will be preserved when plotting them in
## lattice --> all will be first, should be decreasing after
loci.subsets <-
  list(all=function(d)rep(T,nrow(d)))
##        not.all.fixed.same=function(d)!(allcols(d,0)|allcols(d,1)),
##        not.all.fixed=function(d)!allcols(d,c(0,1)),
##        not.all.small.same=function(d)!(allcols(d,lo)|allcols(d,hi)),
##        none.fixed=function(d)onecol(d,c(0,1)),
##        none.small=function(d)onecol(d,lo)&onecol(d,hi)
##        )
read.res.and.subsets <- function(d){
  print(d)
  L <- read.res.tables(d)
  L$subsets <- sapply(paste(d,names(loci.subsets),'/',sep=''),
                      read.res.tables,
                      simplify=F)
  names(L$subsets) <- names(loci.subsets)
  L
}
##debug(read.res.and.subsets)
read.results <- function(dvec)lapply(dvec,read.res.and.subsets)
read.alpha.long <- function(f){
  cutfiles <- c('pop','data','paste')
  out <- paste(f,cutfiles,sep='.')
  names(out) <- cutfiles
  cmd <- paste("cut -b1-5",f,'>',out['pop'],
               '&& cut -b6-',f,'>',out['data'],
               '&& paste',out['pop'],out['data'],'>',out['paste'])
  system(cmd)
  read.table(out['paste'],header=T)
}
##debug(read.alpha.long)
read.nicholson.outfile <- function(f){
  if(!file.exists(f))return(NULL)
  if(length(grep("S.LOC",f)))return(read.s(f))
  if(length(grep("alpha",f)))return(read.alpha.long(f))
  if(length(grep("res_c",f)))return(read.table(f))
  if(length(grep("HMN",f)))return(read.table(f))
  read.table(f,header=T)
}
print.simresult <- function(rl){
  cat("Number of rows in each table:\n")
  print(t(sapply(rl,function(L)sapply(L,dim)[1,])))
}

describe.alleles.s <- function(s){
  types <- unique(s$type)
  m <- sapply(types,function(v)s[s$type==v,])
  res <- sapply(m['s',],
                function(v)summary(factor(v,levels=unique(unlist(m['s',])))))
  colnames(res) <- as.character(types)
  res
}
##debug(describe.alleles.s)
describe.alleles <- function(nsrep){
  s <- results[[rownames(params)[params$N.s.rep==nsrep][1]]]$s
  describe.alleles.s(s)
}
compare.alleles <- function(N){
  ## read data from list structure we read initially
  datatabs <- sapply(names(results[[N]]$subsets),
                     function(l){
                       if(l=='all')L <- results[[N]]$s
                       else L <- results[[N]]$subsets[[l]]$s;
                       L[,-3]},simplify=F)
  ## ghetto filesystem reading method
##   subdir <- paste('results',N,sep='/')
##   sfiles <- dir(subdir,
##                 pattern='S.LOC.txt',
##                 recursive=T,
##                 full.names=T)
##   labs <- gsub('/','',sub("S.LOC.txt",'',sub(subdir,'',sfiles)))
##   labs[labs==''] <- 'all'
##   names(sfiles) <- labs
##   datatabs <- sapply(names(sfiles), function(l)read.s(sfiles[l]), simplify = F)
  svals <- sort(unique(datatabs$all$s))
  typevals <- sort(unique(datatabs$all$type))
  datatabs <- sapply(datatabs,function(L){L$s <- factor(L$s,svals);levels(L$type) <- typevals;L},
                     simplify=F)
  freqtabs <- sapply(names(datatabs),
                     function(l)data.frame(table(datatabs[[l]]),fixed=l),
                     simplify=F)
  proptabs <- lapply(freqtabs,
                     function(df){df$Freq <- df$Freq/freqtabs$all$Freq;df})
  ft <- data.frame(combine(freqtabs),prop=combine(proptabs)$Freq)
  ft$s <- as.numeric(as.character(ft$s))
  ft$S <- paste(ft$type,format(ft$s))
  ##ft$sneg <- ft$s*ifelse(as.character(ft$type)=="BAL",-1,1)
  levs <- sort(unique(ft$S))
  ids <- rev(grep("BAL",levs))
  levs <- levs[c(ids,(ids[1]+1):length(levs))]
  ft$S <- factor(ft$S,levs)
  ft$id <- N
  ft
}
##debug(compare.alleles)
rename <- function(x,m){ # x is names, m is the map
  names(x) <- x
  x[names(m)] <- m
  x
}
