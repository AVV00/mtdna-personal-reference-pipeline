rm(list=ls())
samplename<-grep('.fa$',dir(pattern='.fa'),value=T)
chrM<-read.table(file=samplename,sep='\t',header=F,check.names=F,stringsAsFactors=F)
ref.16569<-toupper(unlist(strsplit(apply(as.matrix(chrM[-1,]),2,function(x){paste(x,collapse="")}),'')))
ref.17119<-c(ref.16569,ref.16569[1:550],rep('',31))
matrix<-matrix(apply(matrix(ref.17119,ncol=70,byrow=T),1,function(x){paste(x,collapse='')}),ncol=1)
res<-rbind('>chrM',matrix)
write.table(res,file=paste0('modified.',samplename),col.names=F,row.names=F,quote=F)

