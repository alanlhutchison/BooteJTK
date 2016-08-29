#!/usr/bin/env Rscript

source("https://bioconductor.org/biocLite.R")
list.of.packages.bioc <- c("limma")
new.packages.bioc <- list.of.packages.bioc[!(list.of.packages.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages.bioc)>0) biocLite(new.packages.bioc,ask=FALSE)

list.of.packages.reg <- c("reshape","reshape2")
new.packages.reg <- list.of.packages.reg[!(list.of.packages.reg %in% installed.packages()[,"Package"])]
if(length(new.packages.reg)>0) install.packages(new.packages.reg,repos='http://cran.us.r-project.org')


# library("optparse")
# option_list = list(
#   make_option(c("-f", "--file"), type="character", default=NULL, 
#               help="dataset file name", metavar="character"),
#   make_option(c("-o", "--out"), type="character", default="out.txt", 
#               help="output file name [default= %default]", metavar="character")
#   make_option(c("-p", "--period"), type="integer", default=24,
#               help="period length [default= %default]", metavar="integer")
# ); 
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# if (is.null(opt$file)){
#   print_help(opt_parser)
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# }
# df = read.table(opt$file, header=TRUE)
# prefix = opt$out
# period = opt$period
library('limma')
library('reshape2')
library('reshape')
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

#fn = 'Gaussian1000_L24.txt'
#pre = 'prefix'
df = read.csv(fn,header=TRUE,sep='\t')
#period = 24
rownames = df[,1]

row.names(df) = rownames
df = df[,-1]
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
colnames(df) <- tx
series = df
MAX = max(table(as.numeric(colnames(series))%%24))
print(MAX)
series.new = NULL
times = c()
t2 = unique(t1%%24)
t2 = t2[ordered(t2)]
print(head(series))
for (h in t2){
  times = c(times,rep(h,dim(series)[1]))
}
for (h in t2){
  ser = series[,as.numeric(colnames(series))%%24==h]
  print(head(ser))
  
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

series.act =  series.new
rownames(series.act) = NULL


series.new['Mean'] = apply(series.act,1,mean)
series.new['SD'] = sqrt(series.ebayes$s2.post)
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

means_out = paste0(pre,'_Means_postLimma.txt')
sds_out   = paste0(pre,'_Sds_postLimma.txt')
ns_out    = paste0(pre,'_Ns_postLimma.txt')

write.table(means,file=means_out,sep='\t',row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(sds,  file=sds_out,sep='\t',  row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(ns,   file=ns_out,sep='\t',   row.names = FALSE, col.names = TRUE,quote=FALSE)

