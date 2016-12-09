#nk.doica=doICA(nk.ica, path , 3)
require(dplyr)
path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/"

nk.ica=matrix_del_rows_scaled(nk,nk.corr,path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")

nk.ica.batch <- doICABatch(nk.ica, seq(2,25,1), path) 
nk.doica <- doICA(nk.ica, path,n=4)

 
b_cells.ica=matrix_del_rows_scaled(b_cells,c(b_cells.cor[[2]],b_cells.cor.2),path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")


caf.ica=matrix_del_rows_scaled(caf,c(),path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")
caf.ica.2=matrix_del_rows_scaled(caf,caf.cor,path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")

t_cells.ica=matrix_del_rows_scaled(t_cells,c(t_cells.cor[[2]],t_cells.cor.2),path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")

macrophages.ica=matrix_del_rows_scaled(macrophages,macrophages.cor,path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")


caf.doica <- running_batch_ica(caf.ica,path)
b_cells.doica <- running_batch_ica(b_cells.ica,path)
caf.doica <- running_batch_ica(caf.ica,path)
caf.doica.2 <- running_batch_ica(caf.ica.2,path)
t_cells.doica <- running_batch_ica(t_cells.ica,path)

macrophages.doica <- running_batch_ica(macrophages.ica,path)
save.image("~/Documents/CURIE2015/R_code/Deconvolution/qc_tirosh081216.RData")

path_global="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/"

macro.doica_8=doICA(macrophages.ica, path_global , 8, corr_folder="CORRELATION_LOCAL_MAX")      
nk.doica_6=doICA(nk.ica, path_global , 6, corr_folder="CORRELATION_LOCAL_MAX")      
caf.ica.2.doica_9=doICA(caf.ica.2, path_global , 9, corr_folder="CORRELATION_LOCAL_MAX")   
caf.ica.doica_9=doICA(caf.ica, path_global , 9, corr_folder="CORRELATION_LOCAL_MAX")
t_cells.ica.doica_3=doICA(t_cells.ica, path_global , 3)
save.image("~/Documents/CURIE2015/R_code/Deconvolution/qc_tirosh081216.RData")
t_cells.ica.doica_3=doICA(t_cells.ica, path_global , 5, corr_folder="CORRELATION_LOCAL_MAX")
save.image("~/Documents/CURIE2015/R_code/Deconvolution/qc_tirosh091216.RData")
####################################################################
# load the annotation database
source("http://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")

library(org.Hs.eg.db)
# set up your query genes
queryGeneNames <- c('WHRN', 'SANS')
queryGeneNames <- c('WHRN')

# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery)
# subset to get your results
result <- aliasSymbol[which(aliasSymbol[,2] == queryGeneNames),5]
result

#########


# 
# hg.rows<-which(aliasSymbol[,2] %in% rownames.test )
# hg.rows.na<-which(!( rownames.test %in% aliasSymbol[,5]  ))
# hg.rows.na[66]
# head(aliasSymbol)
# rownames.test[5]
# aliasSymbol[5,c(2,5)]
# result <- aliasSymbol[which(aliasSymbol[,2] == rownames.test[5]),2:5]
# hugo_rownames.test <- aliasSymbol[hg.rows, c(2,5)]
# length(rownames.test)
# dim(hugo_rownames.test)
# length(unique(hugo_rownames.test[,1]))
# length(hg.rows)
# head(hugo_rownames.test)
# x <- rownames.test



rownames.test <- rownames(t_cells.ica.doica_3[[2]])
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)

hg.rows.na<-which(!( rownames.test %in% aliasSymbol[,5]  ))
correspondance <- lapply(rownames.test[hg.rows.na], function (x) if( x %in% aliasSymbol[,2] == TRUE) {  x = aliasSymbol[which(aliasSymbol[,2] == x),5] } else {x})
hugo <- unlist(correspondance)


rownames.test <- rownames(t_cells.ica.doica_3[[2]])
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)
rownames.test=sub('.', '-', rownames.test, fixed=TRUE)

hg.rows.na<-which(!( rownames.test %in% aliasSymbol[,5]  ))
correspondance <- lapply(rownames.test[hg.rows.na], function (x) if( x %in% aliasSymbol[,2] == TRUE) {  x = aliasSymbol[which(aliasSymbol[,2] == x),5] } else {x})
hugo <- unlist(correspondance)

length(hugo)
length(HUGO[hg.rows.na])
HUGO <- rownames(t_cells.ica.doica_3[[2]])
HUGO[hg.rows.na] <- hugo
hugo  %in%  HUGO[hg.rows.na]

length(hg.rows.na)
tail(hugo,50)
tail(HUGO[hg.rows.na],30)
test1 <- rownames.test[hg.rows.na]


for (i in 1:length(test1)) {
  #print(test1[i] )
  #print(i)
  #print(test1[i] %in% aliasSymbol[,2])
  if( test1[i] %in% aliasSymbol[,2] == TRUE) {
#print(aliasSymbol[which(aliasSymbol[,2] == test1[i]),5])
    test1[i]=aliasSymbol[which(aliasSymbol[,2] == test1[i]),5][1]


  }
}

length(test1)
test1[906]
aliasSymbol[which(aliasSymbol[,2] == "CENPC1"),5]

rownames.test[hg.rows.na][which(rownames.test[hg.rows.na] == "CENPC1")]


# 
# unlist(correspondance ,recursive=FALSE)
# length(unlist(correspondance ))
# rownames.test[22]
# aliasSymbol[22,5]
# which(rownames.test == aliasSymbol[22,2])
# length(hg.rows.na)
# rownames.test[hg.rows.na][1]
