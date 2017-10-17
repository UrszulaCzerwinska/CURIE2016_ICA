# nk.doica
# caf.doica
# b_cells.doica
# t_cells.doica 
# macrophages.doica

name="nk.ica"
path_global
ncomp=4
corr_folder="CORRELATION"
repair <- function(name, path_global, ncomp, corr_folder="CORRELATION") {
  name2=paste(name,"_", ncomp, sep="")

  path_global= paste(path_global, name2, "/", sep="")
  setwd(path_global)
  
  a.imp = paste(paste("A",paste(name, "_numerical.txt", sep = ""),ncomp, sep="_"),"num",sep=".")
  a.imp.path = paste(path_global,a.imp,sep="")
  A=read.delim(a.imp.path, header=FALSE, sep="\t")
  s.imp = paste(paste("S",paste(name, "_numerical.txt", sep = ""),ncomp, sep="_"),"num",sep=".")
  s.imp.path = paste(path_global,s.imp,sep="")
  S=read.delim(s.imp.path, header=FALSE, sep="\t")
  
  ics = paste("IC",c(1:ncomp),sep="")
  colnames(S) = ics
  colnames(A) = ics
  
  rownames(S) = as.matrix(read.delim(paste(path_global,paste(name, "_ids.txt", sep = ""),sep=""), header=FALSE, sep="\t"))
  rownames(A) = as.matrix(read.delim(paste(path_global,paste(name, "_samples.txt", sep = ""),sep=""), header=FALSE, sep="\t"))
  A=A[,1:(ncol(A)-1)]
  S=S[,1:(ncol(S)-1)]
  
  HUGO=change_to_hugo_official(rownames(S))
  
  
  path1 <- paste("../",corr_folder,"/", sep="")
  dir.create(path1)  
  write.table( cbind(HUGO,S), file = paste(path1, paste(name, ncomp,"_S", "xls", sep="."),sep="") , row.names=FALSE, col.names=TRUE, quote= FALSE, sep="\t")
  write.table(data.frame(GENES=rownames(A), A), file = paste(name, ncomp, "full_samples", "txt", sep=".") , row.names=FALSE, col.names=TRUE, quote= FALSE, sep="\t")
  write.table(data.frame(SAMPLES=rownames(S), cbind(S)), file = paste(name,ncomp, "full_genes", "txt", sep=".") , row.names=FALSE, col.names=TRUE, quote= FALSE, sep="\t")
  
  return(list(A,S, HUGO))
}



path_global="/Users/ulala/Documents/CURIE/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/"


nk.doica.4=repair("nk.ica", path_global, ncomp=4)
caf.doica.5=repair("caf.ica", path_global, ncomp=5)
b_cells.doica.3=repair("b_cells.ica", path_global, ncomp=3)
t_cells.doica.3=repair("t_cells.ica", path_global, ncomp=3)
macrophages.doica.2=repair("macrophages.ica", path_global, ncomp=2)

name="t_cells.ica"
ncomp=3

nk.doica.6=repair("nk.ica", path_global, ncomp=6, "CORRELATION_LOCAL_MAX")
caf.doica.9=repair("caf.ica", path_global, ncomp=9, "CORRELATION_LOCAL_MAX")
b_cells.doica.3=repair("b_cells.ica", path_global, ncomp=3, "CORRELATION_LOCAL_MAX")
t_cells.doica.5=repair("t_cells.ica", path_global, ncomp=5, "CORRELATION_LOCAL_MAX")
macrophages.doica.8=repair("macrophages.ica", path_global, ncomp=8, "CORRELATION_LOCAL_MAX")


str(caf.doica)

head(nk.doica[[2]],50)
S.nk=nk.doica[[2]]
dim(S.nk)
test=apply (S.nk, 2, function(x) x[which(x > 0)] = "N/A")

df=nk.doica
corr_folder="CORRELTION_POS_NEG"
split_ic_pos_neg <- function(df, corr_folder) {
  name <- deparse(substitute(df))
  HUGO <- df[[3]]
  S.df.neg <- df[[2]]
  S.df.pos <- df[[2]]
  for (i in 1:ncol(df[[2]]) ) {
    x <- S.df.neg[, i]
    y <- S.df.pos[, i]
    x[which(x >= 0)] = "N/A" #neg
    y[which(y < 0)] = "N/A" #pos
    S.df.neg[,i]=x
    S.df.pos[,i]=y
  }
  path1 <- paste("../",corr_folder,"/", sep="")
  dir.create(path1)  
  
  write.table( cbind(HUGO,S.df.neg), file = paste(path1, paste(name, "neg","_S", "xls", sep="."),sep="") , row.names=FALSE, col.names=TRUE, quote= FALSE, sep="\t")
  write.table( cbind(HUGO,S.df.pos), file = paste(path1, paste(name, "pos","_S", "xls", sep="."),sep="") , row.names=FALSE, col.names=TRUE, quote= FALSE, sep="\t")
  
return(list(S.df.neg,S.df.pos))
  
}

S_pos_neg_nk=split_ic_pos_neg(nk.doica.4, corr_folder="CORRELTION_POS_NEG")
S_pos_neg_caf=split_ic_pos_neg(caf.doica.5, corr_folder="CORRELTION_POS_NEG")
S_pos_neg_b_cells=split_ic_pos_neg(b_cells.doica.3, corr_folder="CORRELTION_POS_NEG")
S_pos_neg_t_cells=split_ic_pos_neg(t_cells.doica.5, corr_folder="CORRELTION_POS_NEG")
S_pos_neg_macrophages=split_ic_pos_neg(macrophages.doica.2, corr_folder="CORRELTION_POS_NEG")

run_correlation(path_global, corr_path = "CORRELTION_POS_NEG")
nk.doica[[2]][,1][which(nk.doica[[2]][,1] > 0)] = "N/A"
#test[which(nk.doica[[2]][,1] > 0)] = "N/A""N/A
nk.doica[[2]]

str(nk.doica)
nk.doica[[3]]

# generate_rnk_files_for_gsea(list(A,S, HUGO), pos_neg = FALSE, rep="rnk/",path_global, name=paste(name,ncomp,sep="."))
list_1 <-  list("S_pos_neg_nk"=S_pos_neg_nk,"S_pos_neg_caf" = S_pos_neg_caf, "S_pos_neg_b_cells" = S_pos_neg_b_cells, "S_pos_neg_t_cells" = S_pos_neg_t_cells, "S_pos_neg_macrophages" = S_pos_neg_macrophages )
for (i in list(S_pos_neg_nk,S_pos_neg_caf,S_pos_neg_b_cells,S_pos_neg_t_cells,S_pos_neg_macrophages)){
  name=deparse(substitute(i))
  generate_rnk_files_for_gsea(i, pos_neg = TRUE, rep="rnk/", path_global) 
  
  
}
lapply(names(list_1),function(x) generate_rnk_files_for_gsea(list_1[[x]], pos_neg = TRUE, rep="rnk/", path_global, name=x))

vec <- c(1:2, 4:20)
for (i in vec){
  repair(paste("METABRIC_0.9", i,sep="_"), "/Users/ulala/Documents/CURIE/Data/breast_cancer_6datasets_ICA/UC/technical/rob_90per_100", 100, corr_folder="CORRELATION")
}
