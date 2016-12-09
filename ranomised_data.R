# GSE72056_melanoma_single_cell <- read.delim("~/Documents/CURIE2015/Data/single.cell.melanoma/GSE72056_melanoma_single_cell.txt", row.names=NULL)
# 
# 
# GSE72056_melanoma_single_cell.w=data.frame(GSE72056_melanoma_single_cell)
row.names(GSE72056_melanoma_single_cell.w)=GSE72056_melanoma_single_cell.w[,1]
GSE72056_melanoma_single_cell.w2=GSE72056_melanoma_single_cell.w[,-1]
GSE72056_melanoma_single_cell.w.t=data.frame(t(GSE72056_melanoma_single_cell.w2))
#
#malignant(1=no,2=yes,0=unresolved)
#non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)


GSE72056_melanoma_single_cell.w[1:10,1:10]
GSE72056_melanoma_single_cell.w.t[1:10,1:10]

library(dplyr)

tumor_cells<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==2, cell_type==0)
dim(tumor_cells)
tumor_cells[1:10,1:10]

t_cells<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==1)
dim(t_cells)

b_cells<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==2)
dim(b_cells)

macrophages<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==3)
dim(macrophages)

endothelial_cells<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==4)
dim(endothelial_cells)

caf<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==5)
dim(caf)

nk<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==6)
dim(nk)

undef<-GSE72056_melanoma_single_cell.w.t%>%add_rownames%>%filter(malignant==1, cell_type==0)
dim(undef)

###simulation
require(dplyr)
t_f=0.6
tc_f=0.1
bc_f=0.05
m_f=0.02
en_f=0.05
c_f=0.015
nk_f=0.015
u_f=0.15
params=data.frame(t_f, tc_f, bc_f, m_f, en_f, c_f, nk_f,u_f)
dataset_list=list(t_cells, b_cells, macrophages,endothelial_cells, caf, nk, undef)
dataset_params_list=list(list(tumor_cells,t_f), list(t_cells,tc_f), list(b_cells,bc_f), list(macrophages,m_f),list(endothelial_cells,en_f), list(caf,c_f), list(nk,nk_f), list(undef, u_f))


nbcell=1000
nsampl=200
#nsampl=3
#apply(params,1, function(x) sample_n(tumor_cells, x, replace=T) )
colnames(simul_samples)=names(res[[1]])

simul_samples=matrix(, nrow = nsampl, ncol = ncol(undef))
#i=1

for (i in 1:nsampl){
  res=lapply(dataset_params_list, function(x) apply(sample_n(x[[1]], x[[2]]*nbcell, replace=T),2, mean))
  df=data.frame(matrix(unlist(res), nrow=length(res), byrow=T),stringsAsFactors=FALSE)
  simul_samples[i,]=apply(df, 2, mean)
}
  
simul_samples.cl= data.frame(simul_samples)%>% select(-(tumor:cell_type))
simul_samples.t=t(scale(simul_samples.cl, scale=FALSE))
row.names(simul_samples.t)=names(res[[1]])[4:23689]

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/simulated_ICA")

write.table(simul_samples.t, file="simul_samples_200_29.08.16_full.txt", quote=FALSE, sep="\t")
write.table(simul_samples.t, file="simul_samples_200_29.08.16_num.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = FALSE)
write.table(row.names(simul_samples.t), file="simul_samples_200_29.08.16_ids.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = FALSE)
write.table(colnames(simul_samples.t), file="simul_samples_200_29.08.16_samples.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names = FALSE)


################################
drp.rows=c("cy80.Cd45.pos.Pd1.neg.S293.E05.S293.comb","cy53.1.CD45.pos.2.D10.S1006.comb","CY88CD45POS_2_D06_S426_comb","CY89A_CD45_POS_10_E06_S246_comb","cy53.1.CD45.pos.2.B12.S984.comb","CY88CD45POS_7_G04_S268_comb","cy60_1_cd_45_pos_4_E08_S56_comb","cy84_Primary_CD45_pos_E08_S440_comb","cy84_Primary_CD45_pos_C07_S415_comb","cy60_1_cd_45_pos_4_C07_S31_comb","cy60_1_cd_45_pos_4_G03_S75_comb","cy84_Primary_CD45_pos_G03_S459_comb","cy60_1_cd_45_pos_4_C11_S35_comb","cy84_Primary_CD45_pos_C11_S419_comb")
NK.sinCell.filtr1<- nk[!(nk$rowname %in% drp.rows),]

###########
drop.rows1=c("CY88CD45POS_2_G03_S459_comb","cy53.1.CD45.pos.2.D10.S1006.comb","CY88CD45POS_2_D06_S426_comb","CY89A_CD45_POS_10_E06_S246_comb","cy53.1.CD45.pos.2.B12.S984.comb","CY88CD45POS_7_G04_S268_comb")
drop.rows2=c("cy60_1_cd_45_pos_4_C08_S32_comb","cy84_Primary_CD45_pos_C08_S416_comb")
drop.rows3=c("cy84_Primary_CD45_pos_E01_S433_comb","cy60_1_cd_45_pos_4_E01_S49_comb")
drp=c(drop.rows1, drop.rows2, drop.rows3)
Macrophages.sinCell.filtr=macrophages[!(macrophages$rowname %in% drp),]
dim(Macrophages.sinCell.filtr)
rm(drop.rows1, drop.rows3, drop.rows2,drp, drp, drp.rows, drp.rows1)
#################################
#rm(list=setdiff(ls(), ls("GSE72056_melanoma_single_cell.w","GSE72056_melanoma_single_cell.w.t","Macrophages.dc_2", "NK.dc")))


###
dim(nk)

mymatrix=matrix(nrow=51, ncol= 100)

bigmatrix=big.matrix(51, 10^10, type="double")

mylist=nk

###########

setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")

counts_plot(macrophages)
counts_plot(nk)
counts_plot(t_cells, ret=TRUE)
counts_plot = function(i, ret = FALSE) {


title=deparse(substitute(i))

mymatrix=matrix(nrow=dim(i)[1], ncol= 100)
z=0
  for (j in seq(0,3.96, 0.04)){
    z=z+1
    #print(j)
  a=dim(i)[2]
  dat=i[,5:a]
  mymatrix[,z]=apply(dat,1,function(x) length(which((x>j)==TRUE)))
  }
rm(dat)

data.tmp=data.frame(cell=i$rowname,mymatrix)
colnames(data.tmp)=c("cell", as.character(seq(0,3.96, 0.04)))
#library(reshape2)
data.tmp.long=melt(data.tmp,id=c("cell"),variable.name = "Threshold")
data.tmp.long$Threshold=as.numeric(data.tmp.long$Threshold)

pl=ggplot(data.tmp.long, aes(x = Threshold, y = value, color = cell, label = cell))+geom_line()+theme_bw()+theme(legend.position="none")+ ylim(0,max(data.tmp.long$value))
rm(i)
ggsave(paste(title, ".pdf", sep=""),pl)

if (ret == TRUE){
  return(data.tmp)
}else{rm(data.tmp)}
rm(mymatrix, data.tmp.long, pl, i, j )

}
