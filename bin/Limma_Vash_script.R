#!/usr/bin/env Rscript


source("https://bioconductor.org/biocLite.R")
list.of.packages.bioc <- c("limma")
new.packages.bioc <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)>0) biocLite(new.packages.bioc,ask=FALSE)


list.of.packages.reg <- c("reshape","reshape2","devtools",'qvalue')
new.packages.reg <- list.of.packages.reg[!(list.of.packages.reg %in% installed.packages()[,"Package"])]
if(length(new.packages.reg)>0) install.packages(new.packages.reg,repos='http://cran.us.r-project.org')

library('limma')
library('reshape2')
library('reshape')
library('qvalue')
library('devtools')

install_github("stephens999/ashr")
library('ashr')

install_github("mengyin/vashr",build_vignettes=TRUE)
library('vashr')


args = commandArgs(trailingOnly=TRUE)

fn = args[1]
pre = args[2]
period = as.numeric(args[3])
## program...
# print(fn)
# print(class(fn))
# print(pre)
# print(class(pre))
# print(period)
# print(class(period))
#fn = '../Hughes_by2-1_head.txt'
#pre = 'prefix'
df = read.csv(fn,header=TRUE,sep='\t')
#period = 24
rownames = df[,1]

row.names(df) = rownames
df = df[,-1]
t1<-(plyr::mapvalues(colnames(df), from=colnames(df), unique(gsub("X", "", colnames(df)))))
t1<-(plyr::mapvalues(t1, from=t1, unique(gsub("ZT", "", t1))))
t1<-as.numeric(t1)
t1<-round((t1-floor(t1))*10*24+floor(t1))
#t1 <- t1 %% period
colnames(df) <- t1
series = df
MAX = max(table(as.numeric(colnames(series))%%24))
print(MAX)
series.new = NULL
times = c()
t2 = unique(t1%%24)
t2 = t2[ordered(t2)]
for (h in t2){
  times = c(times,rep(h,dim(series)[1]))
}
for (h in t2){
  ser = series[,as.numeric(colnames(series))%%24==h]

  
  if (class(ser)=='numeric'){
    for (i in c(1:(MAX-1))){
      ser=cbind(ser,rep(NaN,length(ser)))
    }
  }
    else if(dim(ser)[2]!=MAX){
    for (i in c(1:(MAX-dim(ser)[2]))){
      ser=cbind(ser,rep(NaN,dim(ser)[1]))   

    }
  }
  if (is.null(series.new)){
    series.new = ser
    series.new['ID'] = row.names(series.new)
    row.names(series.new) <- NULL
  }else
  {    
    ser['ID'] = row.names(ser)    
    names(ser) <- names(series.new)
    series.new = rbind(series.new,ser)
  }
}
rownames = series.new$ID
series.new = series.new[,-3]

#rownames =   rownames(series.new)
#rownames(series.new) = NULL
series.new = data.frame(series.new)
print(head(series.new))
print(tail(series.new))

series.fit = limma::lmFit(series.new)
series.ebayes = eBayes(series.fit,robust = TRUE,trend = TRUE)



sehat <- series.ebayes$stdev.unscaled*series.ebayes$sigma
fit.vash <- vash(sehat=sehat, df=series.ebayes$df.residual[1])


series.act =  series.new
rownames(series.act) = NULL


series.new['Mean'] = apply(series.act,1,mean)
#series.new['SD'] = sqrt(series.ebayes$s2.post)
series.new['SD'] = fit.vash$sd.post

series.new['Time'] = times
series.new['ID'] = rownames # gsub('at.*','at',rownames)
series.new['N'] = apply(series.act,1,length)


print('Series.new')
print(head(series.new))
# 
#series.means = series.new[,c('ID','Time','Mean')]
# means =dcast(series.new,ID ~ Time,value.var='Mean')
# sds = dcast(series.new,ID ~ Time,value.var = 'SD')
# ns = dcast(series.new,ID ~ Time,value.var = 'N')

#series.melt = melt(series.new,id.vars=c('ID','Time'))
# means = dcast(series.melt,ID ~ Time,value.var = Mean)
# sds = dcast(series.new,ID ~ Time,value.var = 'SD')
# ns = dcast(series.new,ID ~ Time,value.var = 'N')
series.melt = melt(series.new[,c('ID','Time','Mean','SD','N')],id.vars=c('ID','Time'))

means = dcast(series.melt[series.melt$variable=='Mean',],ID ~ Time ,value.var= 'value')
sds = dcast(series.melt[series.melt$variable=='SD',],ID ~ Time ,value.var= 'value')
ns = dcast(series.melt[series.melt$variable=='N',],ID ~ Time ,value.var= 'value')

print('Means')
print(head(means))

suf='postVash.txt'
means_out = paste0(pre,'_Means_',suf)
sds_out   = paste0(pre,'_Sds_',suf)
ns_out    = paste0(pre,'_Ns_',suf)

write.table(means,file=means_out,sep='\t',row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(sds,  file=sds_out,sep='\t',  row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(ns,   file=ns_out,sep='\t',   row.names = FALSE, col.names = TRUE,quote=FALSE)

