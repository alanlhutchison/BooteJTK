#!/usr/bin/env Rscript


#source("https://bioconductor.org/biocLite.R")
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

#install_github("stephens999/ashr",build_vignettes=FALSE)
library('ashr')

#install_github("mengyin/vashr",build_vignettes=FALSE)
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
#fn = '/Users/alanlhutchison/Desktop/real_data_large/Mauvoisin2/Mauvoisin_prot_top10_jtkready.txt'
#pre = 'prefix'
df = read.csv(fn,header=TRUE,sep='\t')
#period = 24
rownames = df[,1]

row.names(df) = rownames
df = df[,-1]
print(colnames(df))
dfa = read.csv(fn,sep='\t',header=FALSE)
tx = as.vector(t(dfa[1,-1]))
tx<-(plyr::mapvalues(tx, from=tx, gsub("X", "", tx)))
tx<-(plyr::mapvalues(tx, from=tx, gsub("ZT", "", tx)))
tx<-(plyr::mapvalues(tx, from=tx, gsub("CT", "", tx)))
seen= c()
new = c()
for (x in tx)
{
  if (!(x %in% seen)){
    seen <- append(seen,x)
    new <- append(new,x)
  }
  else if (x %in% seen){
    y = x
    while (y %in% tx){
      y <- as.character(as.numeric(y)+24)
      while (y %in% new){
        y <- as.character(as.numeric(y)+24)
      }
    }
    seen<-append(seen,y)
    new <- append(new,y)
  }
}
#t1<-round((t1-floor(t1))*10*24+floor(t1))
#t1 <- t1 %% period
print(tx)
colnames(df) <- as.numeric(tx)
series = df
print(colnames(series))
MAX = max(table(as.numeric(colnames(series))%%24))
print(MAX)
series.new = NULL
times = c()
t2 = unique(as.numeric(tx)%%24)
t2 = t2[ordered(t2)]
for (h in t2){
  times = c(times,rep(h,dim(series)[1]))
}
for (h in t2){
  ser = series[,as.numeric(colnames(series))%%24==h]
  ser.rownames = row.names(series)
  if (class(ser)=='numeric'){
    for (i in c(1:(MAX-1))){
      ser=cbind(ser,rep(NaN,length(ser)))
      #print(ser)
      ser = data.frame(ser)
    }
  }
    else if(dim(ser)[2]!=MAX){
      for (i in c(1:(MAX-dim(ser)[2]))){
        ser=cbind(ser,rep(NaN,dim(ser)[1]))   
      }
      ser = data.frame(ser)
    }
  # else{
  #     print('Leakage')      
  #     print('Class of ser')
  #     print(class(ser))
  #     print('dim(ser)')
  #     print(dim(ser))
  #     print(MAX)
  # 
  # }
  if (is.null(series.new)){
    series.new = ser
    series.new['ID'] = ser.rownames
    row.names(series.new) <- NULL
  }else
  { 
    ser['ID'] = ser.rownames
    names(ser) <- names(series.new)
    series.new = rbind(series.new,ser)
  }
  print(dim(series.new))
}
rownames = series.new$ID
series.new = series.new[,-3]

#print(dim(series.new))
#print(length(times))

#rownames =   rownames(series.new)
#rownames(series.new) = NULL
series.new = data.frame(series.new)
print(head(series.new))
print(tail(series.new))
print('Here we are')
#print(series.new)
#print(class(series.new))
#print(dim(series.new))
#print(dim(as.matrix(series.new,dim(series.new)[1],dim(series.new)[2])))
series.new = data.matrix(series.new)#,dim(series.new)[1],dim(series.new)[2])
#print(class(series.new))
#print(class(series.new[1,1]))
series.fit = limma::lmFit(series.new,na.rm=TRUE)
print('Just did the fit')
#print(series.fit)
for (name in names(series.fit)){
  #print(name)
  if (sum(is.na(series.fit[[name]]))>0){
    #print(name)
    #print(series.fit[[name]])
    series.fit[[name]]<-replace(x<-series.fit[[name]],is.na(x),mean(series.fit[[name]],na.rm=TRUE))
  }
  else{
    series.fit[[name]]<-series.fit[[name]]
  }
}

#series.ebayes = eBayes(series.fit,robust = TRUE,trend = TRUE)



sehat <- series.fit$stdev.unscaled*series.fit$sigma

sehat <- replace(x<-sehat,is.na(x),mean(sehat,na.rm=TRUE))


fit.vash <- vash(sehat=sehat, df=mean(series.fit$df.residual,na.rm=TRUE))

series.new = data.frame(series.new)
series.act =  series.new
rownames(series.act) = NULL
print('series.act')
#print(series.act)
print('mean')
fmean <- function(x){mean(x,na.rm = TRUE)}
series.new['Mean'] = apply(series.act,1,fmean)
#series.new['SD'] = sqrt(series.ebayes$s2.post)
print('sd')
series.new['SD'] = fit.vash$sd.post
print('time')
#print(times)
#print(length(times))
#print(dim(series.new))
series.new['Time'] = times
#print('id')
series.new['ID'] = rownames # gsub('at.*','at',rownames)
#print('n')
series.new['N'] = apply(series.act,1,length)


print('Series.new')
#print(series.new)
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

