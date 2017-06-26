createDF <- function(
  inDirName = '/home/seonwhee/Bioinformatics/IRCR_HTS/'
  
  
){
  
  inFileL = dir(inDirName)
  inFileAll = inFileL[grep('^Nonlin fit of.*.csv',inFileL)]
  inFileAUC = inFileL[grep('^AUC of.*.csv',inFileL)]
  integratedF1 = c()
  integratedF2 = c()
  integratedF3 = c()
  integratedF4 = c()
  
  for (inFN in inFileAll){
    sampN = gsub('_results.csv','',gsub('Nonlin fit of ','',inFN))
    inF_AUCN =inFileAUC[grep(sprintf("AUC of %s_results.csv",gsub("\\+","\\\\+",gsub("\\(","\\\\(",sampN,perl=T),perl=T)),inFileAUC)]
    inF = read.csv(sprintf('%s/%s',inDirName,inFN),sep=",",header=F ,check.names=F)
    inF_AUC = read.csv(file=sprintf('%s/%s',inDirName,inF_AUCN),sep=",",header=F ,check.names=F)
    inF<-t(inF)
    inF_AUC<-t(inF_AUC)
    x=as.numeric(inF_AUC[,3])
    #x[x > 451.5] = NA
    inF_AUC[,3]=x
    
    
    fileLast<- merge(as.matrix(inF[,c(1,2,4,5,6,8)]),as.matrix(inF_AUC[,c(1,3)]),by="V1")
    colnames(fileLast)<-as.vector(fileLast[1,])
    fileLast[(as.numeric(as.character(fileLast[,6])) >1 | is.na(as.numeric(as.character(fileLast[,6])))) | (as.character(fileLast[,2])!="" & as.character(fileLast[,2])!="Hit constraint"),c(3:7)][-1,]<-NA
    
    #HillSlope-integratedF1
    df_HillSlope<-fileLast[-1,c(1,4)]
    colnames(df_HillSlope)= c('Drug',sampN)
    
    #IC50-integratedF2
    df_ic50<-fileLast[-1,c(1,5)]
    colnames(df_ic50)= c('Drug',sampN)
    
    #AUC-integratedF3
    df_AUC<-fileLast[-1,c(1,7)]
    colnames(df_AUC)= c('Drug',sampN)
    
    #logIC50-integratedF4
    df_logIC50<-fileLast[-1,c(1,3)]
    colnames(df_logIC50)= c('Drug',sampN)
    
   
    if (length(integratedF1)==0){
      integratedF1 = df_HillSlope
    } else {      integratedF1= merge(integratedF1,df_HillSlope,all=T)
    }
    if (length(integratedF2)==0){
      integratedF2 = df_ic50
    } else {
      integratedF2= merge(integratedF2,df_ic50,all=T)
    }
    
    if (length(integratedF3)==0){
      integratedF3 = df_AUC
    } else {
      integratedF3= merge(integratedF3,df_AUC,all=T)
    }
    if (length(integratedF4)==0){
      integratedF4 = df_logIC50
    } else {
      integratedF4= merge(integratedF4,df_logIC50,all=T)
      
    }
  }
  
  
  outFileName1 = '/home/seonwhee/Bioinformatics/IRCR_HTS/HillSlope.csv'
  outFileName2 = '/home/seonwhee/Bioinformatics/IRCR_HTS/IC50.csv'
  outFileName3 = '/home/seonwhee/Bioinformatics/IRCR_HTS/AUC.csv'
  outFileName4 = '/home/seonwhee/Bioinformatics/IRCR_HTS/logIC50.csv'
  write.table(t(integratedF1),file=outFileName1,sep=",",quote=F, col.names=F,row.names=T)
  write.table(t(integratedF2),file=outFileName2,sep=",",quote=F, col.names=F,row.names=T)
  write.table(t(integratedF3),file=outFileName3,sep=",",quote=F, col.names=F,row.names=T)
  write.table(t(integratedF4),file=outFileName4,sep=",",quote=F, col.names=F,row.names=T)
}

createDF(inDirName='/home/seonwhee/Bioinformatics/IRCR_HTS/') 

