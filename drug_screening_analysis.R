co.var <- function(x,na.rm=TRUE) 100*(sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))

############   Replace Rownames by Row numbers ###########
get_Index <- function(DF, rname){
  rvalue = which(rownames(DF)==rname)
  return(rvalue)
}

############# 행 번호를 알파벳으로 매긴다 #################
data_Preprocess <- function(inFileName){
  inFile = read.csv(inFileName,sep=",",header=T, row.names=NULL, check.names=F,skip=9,nrow=16)
  rownames(inFile) <- paste("", LETTERS[1:16], sep = "") 
  inFile[,1] <- NULL
  return(inFile)
}
############## Fetching Negative Control data ############  
get_DMSO <- function(inputData){
  DMSO <- inputData[get_Index(inputData,"I"):get_Index(inputData, "O"),3]  
  return(DMSO)
}
############## Fetching Positive Control data #############
get_Positive_Control <- function(inputData){
  Pctrl <- inputData[get_Index(inputData,"B"):get_Index(inputData, "G"),2]  
  return(Pctrl)
}
############### Calculating Z-factor ######################
get_Z_factor <- function(muP, muN, sigmaP, sigmaN){
  z_factor <- 1 - 3 * ((sigmaN+sigmaP)/abs(muN-muP))  
  return(z_factor)
}
###########################################################
get_Screening_data <- function(inFile, DF){
  dmso = get_DMSO(inFile)
  avg.dmso = mean(dmso,na.rm=T)
  norm.DF = DF/avg.dmso*100
  return(norm.DF)
}
###########################################################
get_Control_stats <- function(inFile){
  dmso = get_DMSO(inFile)
  cv = co.var(dmso)
  sd.dmso = sd(dmso,na.rm=T) 
  avg.dmso = mean(dmso,na.rm=T)
  pos_control = get_Positive_Control(inFile)
  sd.pos_control = sd(pos_control, na.rm = T)
  avg.pos_control = mean(pos_control, na.rm = T)
  Z_factor = get_Z_factor(avg.pos_control, avg.dmso, sd.pos_control, sd.dmso)
  results = c(cv, avg.dmso, Z_factor)
  return(results)
}
############################################################
matching_Drug_by_plate <- function(normDF1, normDF2, PANEL, rowFrom, rowTo){
  outDirName = '/home/seonwhee/Bioinformatics/IRCR_HTS/'
  reoutDirName = outDirName
  for (i in length(PANEL)){
    pair1 <- normDF1[get_Index(normDF1, rowFrom):get_Index(normDF1, rowTo),i] 
    pair2 <- normDF2[get_Index(normDF2, rowFrom):get_Index(normDF2, rowTo),i]
    drug.df <- data.frame(pair1, pair2, row.names = NULL)
    drugName <- PANEL[i]
    colnames(drug.df) <- c(drugName, drugName)
    outFileName = sprintf('%s/%s.csv',outDirName,drugName)
    if (file.exists(outFileName)){
      out.df = read.csv(outFileName,header=T,check.names=F)
      out.df = cbind(out.df, drug.df)
    } else{
      out.df = drug.df
    }
    write.table(out.df,file=outFileName,sep=",",quote=F, col.names=T,row.names=FALSE)
  }
}
#############################################################
GBM <- function(inFile1,  inFile2, plate) 
{
    if (plate == 1) {
      df1 = inFile1[-c(1,16),-c(1:3,24:26)]
      df2 = inFile2[-c(1,16),-c(1:3,24:26)]
      panel_A <- c("AZD2014", "AZD3759", "AZD4547", "AZD5363", "AZD6094", "AZD6738", "AZD9291",
                   "Rapamycin (Sirolimus)", "AZD1775 (MK-1775)", "Simvastatin", "triptolide", "Temozolomide",
                   "Imatinib Mesylate (STI571)", "Pazopanib (GW786034)", "Sunitinib Malate", 
                   "Olaparib (AZD2281, KU0059436)", "Regorafenib (BAY 73-4506)", "Erlotinib", "Gefitinib",
                   "Afatinib")
      panel_B <- c("BEZ235", "BGJ398", "BKM120 (NVP-BKM120)", "BYL719", "Cediranib (AZD2171)", "Ceritinib (LDK378)",
                   "Crizotinib (PF-02341066)", "Dacomitinib", "Everolimus (RAD001)", "Ibrutinib (PCI-32765)", 
                   "Lapatinib", "LY2835219", "Neratinib (HKI-272)", "Panobinostat", "PD0332991 (Palbociclib) HCl",
                   "Selumetinib (AZD6244)", "Sorafenib", "Trametinib", "Vandetanib (ZD6474)", "Vemurafenib")  }
    else if(plate == 2){
      df1 = inFile1[-c(1,16),-c(1:3,19:26)]
      df2 = inFile2[-c(1,16),-c(1:3,19:26)]
      panel_A <- c("ABT-199 GDC-0199", "ABT-888 (Veliparib, NSC737664)", "alectinib", "Avagacestat (BMS-708163)",
                    "Axitinib (AG-013736)", "Bosutinib (SKI-606)", "Cabozantinib (XL184)", "Dabrafenib",
                    "Dasatinib (BMS-354825)", "decitabine", "defactinib", "Dovitinib (TKI258, CHIR258)",
                    "Fluorouracil (5-Fluoracil, 5-FU)", "Foretinib (XL880)", "Galunisertib (LY2157299)")
      panel_B <- c("Idelalisib", "Irinotecan", "LY2109761", "Niclosamide", "Nilotinib (AMN-107)", "NVP-AEW541",
                    "Oxaliplatin", "PF-05212384 (PKI-587)", "PLX3397", "Quizartinib (AC220)", 
                    "Saracatinib (AZD0530)", "Semagacestat (LY450139)", "TGX-221", "tozasertib", "XAV939")
    }
  norm.df1 = get_Screening_data(inFile1, df1)
  norm.df2 = get_Screening_data(inFile2, df2)
  matching_Drug_by_plate(norm.df1, norm.df2, panel_A, "B", "H")
  matching_Drug_by_plate(norm.df1, norm.df2, panel_B, "I", "O")
}
####################################################################################################################
################################# Code Execution Part ###################################################
cellName = 'MBT17_341T_M18'
qcFileName = '/home/seonwhee/Bioinformatics/IRCR_HTS/QC.txt'
inFile1 = data_Preprocess('/home/seonwhee/Bioinformatics/IRCR_HTS/3.csv')
inFile2 = data_Preprocess('/home/seonwhee/Bioinformatics/IRCR_HTS/4.csv')
inFile3 = data_Preprocess('/home/seonwhee/Bioinformatics/IRCR_HTS/1.csv')
inFile4 = data_Preprocess('/home/seonwhee/Bioinformatics/IRCR_HTS/2.csv')

GBM(inFile1, inFile2, 1)
GBM(inFile3, inFile4, 2)
qcDF = data.frame(get_Control_stats(inFile1), get_Control_stats(inFile2), get_Control_stats(inFile3),
                  get_Control_stats(inFile4))
rownames(qcDF) <- c("CV", "DMSO(AVE)", "Z-Factor")
colnames(qcDF) <- c("plate 1", "plate 2", "plate 3", "plate 4")
View(qcDF)
write.table(qcDF, file=qcFileName, append=F, sep = "\t")