selectHit <- function(inFileName, outDirName, cutoff, drugListName){
  library(ggplot2)
  library(gridExtra)
  require(grid)
  library(gtable)
  library(dplyr)
  library(stringr)
  
  ### qcdf = read.csv('/home/seonwhee/Bioinformatics/IRCR_HTS/QC.txt', sep = "\t", header = T) # QC data.frame OPTIONAL!!!!
  
  ### Preprocess the overall AUC data table
  set.seed(1)
  df1 = read.csv(inFileName,sep=",",header=T,check.names=F) # read AUC.csv
  df = t(df1)                                               # Transpose --> columns : sample names, Row indices : drug names
  colnames(df)=as.character(df[1,])                         # 1st row will be a header of the table
  df=df[-1,]                                                # delete the 1st row which became a header
  df = data.frame("Drug"=rownames(df),df,check.names=F)     # add a column named "Drug"
  write.table(df,file=sprintf('%soriginalDF.txt',outDirName),sep="\t",quote=F, col.names=T,row.names=T)
  dlm = read.delim(drugListName,check.names=F)              # read a delimited text file(170406_druglist_plot.txt)
  # auc_df contains colnames : "Drug" "Drug2" "Target" "Color" "BT17_020T_M8_SNU" "BT17_1342T_M8" .....
  auc_df = merge(dlm,df,by.x=names(dlm)[1],by.y=names(df)[1],sort=F) # merge two data.frames that share the common rownames
  #className = c("red", "orange", "purple", "green", "yellow", "grey", "pink", "black", "white")
  dat = c()

  for (cell in colnames(auc_df)[-(1:4)]){
    dat = rbind(dat,cbind(auc_df[,1:4],"auc"=as.numeric(as.character(auc_df[,cell])),"Cell"=cell, "label_colors"="black"))
    # create a list that has colnames: "Drug" "Drug2" "Target" "Color" "auc" "Cell" "label_colors"
  }
  ############### Cutoff several data ####################################
  dat = dat[dat$auc < 451.5,]                                 # only get AUC values less than 451.5
  dat = dat[!is.na(dat$auc),]                                 # remove NA
  dat = as.data.frame(dat)
  ####### 위 과정에서 sampleID별 dimension이 달라진다 Boxplot그리기 위해서는 data를 다시 정리해야 한다##################
  
  dat <- within(dat, {
    Drug2 = factor(dat$Drug2,levels=(dlm$Drug2),ordered=T)
    auc = as.numeric(as.character(auc))
    grp = 'zcontrol'
    hit = "no"          # initialize hit select
    label_colors = as.character(dat$label_colors)
    label_colors = "black"
    
  })
  cellL = colnames(df)[-1]    # column "cell"
  ################################################ sampleID별 data subsetting#########################################################################
  for (sampleID in cellL){
    x=as.matrix(auc_df[,-(1:4)])                              # convert 63 by 6 list to 63 by 6 matrix
    mode(x) = "numeric"
    print(sampleID)
    znorm_df = cbind(as.character(auc_df[,1]), as.character(auc_df[,2]), as.character(auc_df[,3]),t(apply(x,1,function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T))))
    hit = which(as.numeric(znorm_df[,sampleID])<= cutoff)      # select rownumbers that filtered by cutoff value
    print(length(hit))
    if (length(hit) == 0){
      hit_drug <- "None"
      hit_value <- "None"
      hit_target <- "No Recommended Drugs"
      hit_df = data.frame("Drug"=hit_drug,"Target"=hit_target,"z_score"= hit_value, stringsAsFactors = F)
    }
    else {
      hit_drug = znorm_df[hit,2]                                # znorm_df[hit,1] --> "Drug",  znorm_df[hit,2] --> "Drug2"  
      hit_value = as.numeric(znorm_df[hit,sampleID])
      hit_target = znorm_df[hit,3]
      hit_df = data.frame("Drug"=hit_drug,"Target"=hit_target,"z_score"= sprintf("%.2f", hit_value), stringsAsFactors = F)
      hit_df = hit_df[order(hit_df$z_score,decreasing=T),]
    }
    
    for(i in 1:length(hit_df$Target)){
      if(str_length(hit_df[i,"Target"]) > 30){
        hit_df[i,"Target"] = str_c(as.character(str_sub(hit_df[i,"Target"], end = 18)), "\n", as.character(str_sub(hit_df[i,"Target"], start = 19, end = -1)))
      }
    }

    subDF <- subset(dat, dat$Cell==sampleID)
    subDF <- within(subDF, {
      auc = as.numeric(as.character(auc))
      grp = sampleID
      hit[Drug2 %in% hit_drug] = "select"
      label_colors[Drug2 %in% hit_drug] = "red"
      
    })
    subDF = merge(x=subDF,y=df,by.x=names(subDF)[1], sort=F, all.x=T) # merge two data.frames that share the common rownames
    
    newList = c()
    for (cell in colnames(subDF)[-(1:9)]){
      newList = rbind(newList, cbind(subDF[,1:9],"AUC"=as.numeric(as.character(subDF[,cell]))))
    }
    
    newDF = as.data.frame(newList)
    
    col = as.character(dlm$Color)  # color discrimination by class
    
    #### Drug Class ################# user-editing part
    levels(newDF$Color)[levels(newDF$Color)=="red"] = "Class1"
    levels(newDF$Color)[levels(newDF$Color)=="orange"] = "Class2"
    levels(newDF$Color)[levels(newDF$Color)=="purple"] = "Class3"
    levels(newDF$Color)[levels(newDF$Color)=="green"] = "Class4"
    levels(newDF$Color)[levels(newDF$Color)=="yellow"] = "Class5"
    levels(newDF$Color)[levels(newDF$Color)=="blue"] = "Class6"
    levels(newDF$Color)[levels(newDF$Color)=="white"] = "Class7"
    levels(newDF$Color)[levels(newDF$Color)=="cyan"] = "Class8"
    levels(newDF$Color)[levels(newDF$Color)=="pink"] = "Class9"
    levels(newDF$Color)[levels(newDF$Color)=="grey"] = "Class10"
    levels(newDF$Color)[levels(newDF$Color)=="black"] = "Class11"
    
    ## code for convinience 
    #ind = 1
    #for(individual_col in unique(col)){
    #  levels(newDF$Color)[levels(newDF$Color)==individual_col] <- paste(c("Class", as.character(ind)), collapse = '')
    #  ind = ind + 1
    #}
    
    #################################BOX PLOT Drawing################################################
    main_title <- textGrob(sampleID, gp=gpar(fontsize=18))
    
    idinfoL <- c("Drug plate ver  :", "Reference set  :")
    
    idinfo_inputL <- c("default", "default") ### user-editing part
    infoL_df <- data.frame(idinfoL, idinfo_inputL)
    infoL_tbl <- tableGrob(infoL_df, rows = NULL, cols = NULL, theme = ttheme_minimal(base_size = 10))
    
    idinfoR <- c("Screening date  :", "Analysis date  :", "QC result  :")
    idinfo_inputR <- c("2017-06-01", "2017-06-01", "Pass: CV<30%, DMSO>10,000") ### user-editing part
    infoR_df <- data.frame(idinfoR, idinfo_inputR)
    infoR_tbl <- tableGrob(infoR_df, rows = NULL, cols = NULL, theme = ttheme_minimal(base_size = 10))
    
    arrange_A <- arrangeGrob(infoL_tbl, infoR_tbl, widths=c(2/5, 3/5), ncol=2)
    arrange_Ratio <- c()
    boxplot_legend = "bottom"
    ################### Hit Selection Table #########################
    tbl_title =  textGrob(sprintf("Recommended Drugs (Z < %.1f)", cutoff), gp=gpar(fontsize=16))
    if (nrow(hit_df) > 20){
      if(nrow(hit_df) > 35){
        tsize = 5
      }
      else{
        tsize = 6
      }
      tbl_sub1 = tableGrob(hit_df[1:19,], rows = NULL, theme = ttheme_default(base_size = tsize))
      tbl_sub2 = tableGrob(hit_df[20:nrow(hit_df),], rows = NULL, theme = ttheme_default(base_size = tsize))
      tbl = arrangeGrob(tbl_sub1, tbl_sub2, widths = c(1/2, 1/2), ncol = 2)
      drugs_tbl = arrangeGrob(tbl_title, tbl, heights = c(1/20, 19/20), ncol = 1)
      arrange_Ratio <- c(2.5/29.7, 2/29.7, 12.5/29.7, 12.7/29.7)
      boxplot_legend = "right"
    }
    else {
      tbl_padding = unit(5, "mm")
      if (nrow(hit_df) < 7){
        tbl = tableGrob(hit_df, rows = NULL, theme = ttheme_default(base_size = 15))
        arrange_Ratio <- c(2.5/29.7, 2/29.7, 16.5/29.7, 8.7/29.7)
      }
      else if (nrow(hit_df) > 16){
        tbl = tableGrob(hit_df, rows = NULL, theme = ttheme_default(base_size = 9))
        arrange_Ratio <- c(2/29.7, 2/29.7, 12/29.7, 13.7/29.7)
        boxplot_legend = "right"
      }
      else {
        tbl = tableGrob(hit_df, rows = NULL, theme = ttheme_default(base_size = 12))
        arrange_Ratio <- c(2.5/29.7, 2/29.7, 12.5/29.7, 12.7/29.7)
      }
      tbl <- gtable_add_rows(tbl, heights = grobHeight(tbl_title)+tbl_padding, pos = 0)
      drugs_tbl <- gtable_add_grob(tbl, tbl_title, 1, 1, 1, ncol(tbl))
    }
    
    ################ Defining Box plot #######################
    bp=ggplot(data=newDF, aes(x=Drug2, y=AUC, fill=factor(Color), colors = factor(Color))) + geom_boxplot() + xlab("") +
      guides(fill=guide_legend(title = "Drug Class"))
    bp = bp + theme(axis.text.x=element_text(angle=90,hjust=1,color=subDF$label_colors)) 
    bp = bp + theme(legend.position = boxplot_legend)
    bp = bp + geom_point(aes(x=Drug2,y=auc,col=hit,shape=grp, size=hit))+
      scale_colour_manual(values=c('black','red'),guide=F)+
      scale_shape_manual(values=c(4,32),guide=F)+
      scale_size_manual(values=c(2,2),guide=F)
    
    ############## Layout and Save jpeg ############################
    Overall_arrange = arrangeGrob(main_title, arrange_A, bp, drugs_tbl, heights = arrange_Ratio, ncol=1)
    ggsave(filename=sprintf('%s/%s_auc_plot.jpeg',outDirName,sampleID),plot=Overall_arrange, width = 21, height = 29.7, units = "cm")
    write.table(hit_df, sprintf("%s/%s_hit_drug.txt",outDirName,sampleID),sep="\t",quote=F, col.names=T,row.names=F)
  }

}

# CODE Execution
selectHit('/home/seonwhee/Bioinformatics/IRCR_HTS/Input_dir/AUC.csv','/home/seonwhee/Bioinformatics/IRCR_HTS/Output_dir',-0.5,'/home/seonwhee/Bioinformatics/IRCR_HTS/170530_druglist_plot.txt')
