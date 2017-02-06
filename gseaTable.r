setwd("/Users/ulala/Documents/CURIE/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/GSEA/")
tab=NULL
count=0
for (i in list.dirs(".", recursive = FALSE)) {
  print(i)
  count=count+1
  print(count)
  require(qpcR)
  for (f in list.files(i, pattern="top_gsea_report_for_na_neg*")) {

      if (file.info(paste(i,f,sep="/"))$size == 0) {
        column_neg <- data.frame("NA")
      }
      else {
        column_neg <- read.table(paste(i,f,sep="/"))
      }
    #print(paste(i,f,sep="/"))
    colnames(column_neg)=paste(i, "NEG", sep="_")
    #print(paste(i, "NEG", sep="_"))
    #print(colnames(column_neg))
  }
  
  for (f in list.files(i, pattern="top_gsea_report_for_na_pos*")) {
    
    if (file.info(paste(i,f,sep="/"))$size == 0) {
      column_pos <- data.frame("NA")
    }
    else {
      column_pos <- read.table(paste(i,f,sep="/"))
    }
    #print(paste(i,f,sep="/"))
    colnames(column_pos)=paste(i, "POS", sep="_")
    #print(paste(i, "POS", sep="_"))
    #print(colnames(column_pos))
  }
  
  pos_neg <- qpcR:::cbind.na(column_neg, column_pos )
  tab <- qpcR:::cbind.na(tab, pos_neg)
  print(pos_neg)
  
}
remove(pos_neg)
write.table(tab[,2:ncol(tab)], file="./top_genes_tr6.txt",quote=FALSE, sep="\t", row.names = FALSE)
#dftest=read.table("./t_cells.ica.3.IC1.rnk_gsea.GseaPreranked.1483384158331/top_gsea_report_for_na_neg_1483384158331.txt")
#colnames(dftest)="tot"
head(  tab)


tab2=NULL
count=0
for (i in list.dirs(".", recursive = FALSE)) {
  print(i)
  count=count+1
  print(count)
  require(qpcR)
  for (f in list.files(i, pattern="top_tr3_gsea_report_for_na_neg*")) {
    
    if (file.info(paste(i,f,sep="/"))$size == 0) {
      column_neg <- data.frame("NA")
    }
    else {
      column_neg <- read.table(paste(i,f,sep="/"))
    }
    #print(paste(i,f,sep="/"))
    colnames(column_neg)=paste(i, "NEG", sep="_")
    #print(paste(i, "NEG", sep="_"))
    #print(colnames(column_neg))
  }
  
  for (f in list.files(i, pattern="top_tr3_gsea_report_for_na_pos*")) {
    
    if (file.info(paste(i,f,sep="/"))$size == 0) {
      column_pos <- data.frame("NA")
    }
    else {
      column_pos <- read.table(paste(i,f,sep="/"))
    }
    #print(paste(i,f,sep="/"))
    colnames(column_pos)=paste(i, "POS", sep="_")
    #print(paste(i, "POS", sep="_"))
    #print(colnames(column_pos))
  }
  
  pos_neg2 <- qpcR:::cbind.na(column_neg, column_pos )
  tab2 <- qpcR:::cbind.na(tab2, pos_neg2)
  print(pos_neg)
  
}
head(tab2)
write.table(tab2[,2:ncol(tab2)], file="./top_genes_tr3.txt",quote=FALSE, sep="\t", row.names = FALSE)
