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

if (length(args)>3){
  bool.rnaseq = TRUE
  }

#fn = '~/Desktop/Allada_raw/adjcounts/by2/reps/HC18_FatBody_ZT_ADJcounts_jtkready.txt'
#pre = 'Allada_HC18_FatBody'
#period = 24
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
colnames(df) <- as.numeric(tx)
series = df
MAX = max(table(as.numeric(colnames(series))%%24))
print(MAX)
series.new = NULL
times = c()

t2 = unique(as.numeric(tx)%%24)
t2 = sort(t2)
print(head(series))
for (h in t2){
  times = c(times,rep(h,dim(series)[1]))
}
for (h in t2){
  if ( is.null(dim(series[,as.numeric(colnames(series))%%24==h])) ){
    ser.voom = voom(ser)
    ser.voom$E = series[,as.numeric(colnames(series))%%24==h]
    ser.voom$weights = rep(mean(ser.voom$weights),dim(ser.voom$weights)[1])
 }else{
  ser = series[,as.numeric(colnames(series))%%24==h]
  print(head(ser))
  if (bool.rnaseq){ser.voom = voom(ser)  }else{ser.voom = vooma(ser)}
  }
  if (class(ser.voom$E)=='numeric'){
    for (i in c(1:(MAX-1))){
      ser.voom$weights = cbind(ser.voom$weights,rep(NaN,length(ser)))
      ser.voom$E = cbind(ser.voom$E,rep(NaN,length(ser)))
      
      #ser=cbind(ser,rep(NaN,length(ser)))
      #ser=data.frame(ser)
    }
  }
    else if(dim(ser.voom$E)[2]!=MAX){
    for (i in c(1:(MAX-dim(ser)[2]))){
      ser.voom$weights = cbind(ser.voom$weights,rep(NaN,dim(ser)[1]))      
      ser.voom$E = cbind(ser.voom$E,rep(NaN,dim(ser)[1]))      
      
      #ser=cbind(ser,rep(NaN,dim(ser)[1]))   
      #ser = data.frame(ser)
    }
  }
  print('Check out ser here!')
  print(head(ser))
  if (is.null(series.new)){
    series.new = ser.voom
    #series.new['ID'] = row.names(series.new)
    rownames = row.names(series.new$weights)
    
    #row.names(series.new) <- NULL
  }else
  {    
    #ser['ID'] = row.names(ser)    
    rownames = c(rownames,row.names(series.new$weights))
    #series.new = rbind(series.new,ser)
    names(ser.voom$weights) <- names(series.new$weights)    
    series.new$weights = rbind(series.new$weights,ser.voom$weights)
    names(ser.voom$E) <- names(series.new$E)    
    series.new$E = rbind(series.new$E,ser.voom$E)
    
  }
}
#rownames = series.new$ID
#series.new = series.new[,-3]
series.new$design  = NULL
series.new$targets = NULL
#rownames =   rownames(series.new)
#rownames(series.new) = NULL
#series.new = data.frame(series.new)
#print(head(series.new))
#print(tail(series.new))

#remove.packages('limma')
#biocLite('limma')
#library(limma)
####### HERE WE DO THE VOOM THING #######
### ALREADY DID IT ####

####### HERE WE DO THE EBAYES THING ########
series.fit = limma::lmFit(series.new)
series.ebayes = limma::eBayes(series.fit,robust = TRUE,trend = TRUE)

series.act =  series.new
rownames(series.act) = NULL

series.new <- data.frame(Mean=rep(1.1,dim(series.act$E)[1]),
                 SD=rep(1.1,dim(series.act$E)[1]), 
                 SDpre=rep(1.1,dim(series.act$E)[1]), 
                 Time = rep(1.1,dim(series.act$E)[1]),
                 ID = rep('A',dim(series.act$E)[1]),
                 N = rep(1.1,dim(series.act$E)[1]),
                 stringsAsFactors=FALSE) 
fmean <- function(x){mean(x,na.rm = TRUE)}
series.new['Mean'] = apply(series.act$E,1,fmean)
series.new['SD'] = sqrt(series.ebayes$s2.post)
series.new['SDpre'] = series.ebayes$sigma
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
series.melt = melt(series.new[,c('ID','Time','Mean','SD','SDpre','N')],id.vars=c('ID','Time'))

means = dcast(series.melt[series.melt$variable=='Mean',],ID ~ Time ,value.var= 'value')
sds = dcast(series.melt[series.melt$variable=='SD',],ID ~ Time ,value.var= 'value')
ns = dcast(series.melt[series.melt$variable=='N',],ID ~ Time ,value.var= 'value')
sdspre = dcast(series.melt[series.melt$variable=='SDpre',],ID ~ Time ,value.var= 'value')

print('Means')
print(head(means))

means_out = paste0(pre,'_Means_postLimma.txt')
sds_out   = paste0(pre,'_Sds_postLimma.txt')
ns_out    = paste0(pre,'_Ns_postLimma.txt')
sdspre_out = paste0(pre,'_Sds-pre_postLimma.txt')


write.table(means,file=means_out,sep='\t',row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(sds,  file=sds_out,sep='\t',  row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(ns,   file=ns_out,sep='\t',   row.names = FALSE, col.names = TRUE,quote=FALSE)
write.table(sdspre,  file=sdspre_out,sep='\t',  row.names = FALSE, col.names = TRUE,quote=FALSE)



