---
  title: "QC Tirosh et al."
output: html_notebook
---
  
  ```{r loading data}
load("~/Documents/CURIE2015/R_code/miunimal_melanoma.RData")

```

```{r data description}
# GSE72056_melanoma_single_cell.w=data.frame(GSE72056_melanoma_single_cell)
#row.names(GSE72056_melanoma_single_cell.w)=GSE72056_melanoma_single_cell.w[,1]
#GSE72056_melanoma_single_cell.w2=GSE72056_melanoma_single_cell.w[,-1]
#GSE72056_melanoma_single_cell.w.t=data.frame(t(GSE72056_melanoma_single_cell.w2))
#
#malignant(1=no,2=yes,0=unresolved)
#non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)
GSE72056_melanoma_single_cell.w[1:10,1:10] #samples in columns
GSE72056_melanoma_single_cell.w.t[1:10,1:10] #samples in rows
```

```{r data split into subsets (labelled)}
require(dplyr)

tumor_cells <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  2, cell_type == 0)

dim(tumor_cells)
tumor_cells[1:10, 1:10]

t_cells <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 1)
dim(t_cells)

b_cells <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 2)

dim(b_cells)

macrophages <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 3)

dim(macrophages)

endothelial_cells <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 4)

dim(endothelial_cells)

caf <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 5)

dim(caf)

nk <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 6)

dim(nk)

undef <-
  GSE72056_melanoma_single_cell.w.t %>% add_rownames %>% filter(malignant ==
                                                                  1, cell_type == 0)

dim(undef)

```
```{r counts plot function}

counts_plot <- function(i, ret = FALSE) {
  title <- deparse(substitute(i))
  
  mymatrix <- matrix(nrow = dim(i)[1], ncol = 100)
  z <- 0
  for (j in seq(0, 3.96, 0.04)) {
    z <- z + 1
    #print(j)
    a <- dim(i)[2]
    dat <- i[, 5:a]
    mymatrix[, z] <- apply(dat, 1, function(x)
      length(which((x > j) == TRUE)))
  }
  rm(dat)
  
  data.tmp = data.frame(cell = i$rowname, mymatrix)
  colnames(data.tmp) = c("cell", as.character(seq(0, 3.96, 0.04)))
  require(reshape2)
  data.tmp.long = melt(data.tmp, id = c("cell"), variable.name = "Threshold")
  data.tmp.long$Threshold = as.numeric(data.tmp.long$Threshold)
  
  require(ggplot2)
  pl <- ggplot(data.tmp.long,
               aes(
                 x = Threshold,
                 y = value,
                 color = cell,
                 label = cell
               )) + geom_line() + theme_bw() + theme(legend.position = "none") + ylim(0, max(data.tmp.long$value))
  
  rm(i)
  
  ggsave(paste(title, ".pdf", sep = ""), pl)
  
  if (ret == TRUE) {
    return(list(data.tmp,data.tmp.long))
  } else{
    rm(data.tmp)
  }
  
  #rm(mymatrix, data.tmp.long, pl, i, j)
  
}

```

```{r testing plotting and cells deletion}
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")
nk_plot <- counts_plot(nk, ret=TRUE)


library(dplyr)
delete.top    <- nk[[2]] %>% filter(Threshold == '50') %>% arrange(value)%>%head(3)%>% select(cell)
delete.bottom <- nk[[2]] %>% filter(Threshold == '50') %>% arrange(value)%>%tail(4)%>% select(cell)
nk_curated<-subset(nk[[2]], !(cell %in% as.matrix(delete.top) | cell %in% as.matrix(delete.bottom) ))
pl <- ggplot(nk_curated,
             aes(
               x = Threshold,
               y = value,
               color = cell,
               label = cell
             )) + geom_line() + theme_bw() + theme(legend.position = "none") + ylim(0, max(nk_curated$value))

nk_plot[[2]][,1:3]

df <- data.frame(x = c(10, 4, 1, 6, 3, 1, 1))
df %>% top_n(2)

# Negative values select bottom from group. Note that we get more
# than 2 values here because there's a tie: top_n() either takes
# all rows with a value, or none.
df %>% top_n(-2)



```

```{r delete cells}





delete_cells <-
  function(df.long,
           thr = 30,
           top = 3,
           bottom = 3,
           path,
           long = FALSE) {
    title <- deparse(substitute(df.long))
    df.long = data.frame(df.long)
    delete.top    <-
      df.long %>% filter(Threshold == thr) %>% arrange(value) %>% head(bottom) %>% select(cell)
    delete.bottom <-
      df.long %>% filter(Threshold == thr) %>% arrange(value) %>% tail(top) %>% select(cell)
    curated <-
      subset(df.long, !(
        cell %in% as.matrix(delete.top) |
          cell %in% as.matrix(delete.bottom)
      ))
    require(ggplot2)
    
    pl <- ggplot(curated,
                 aes(
                   x = Threshold,
                   y = value,
                   color = cell,
                   label = cell
                 )) + geom_line() + theme_bw() + theme(legend.position = "none") + ylim(0, max(curated$value))
    
    setwd(path)
    
    ggsave(paste(title, "curated.pdf", sep = "_"), pl)
    
    cells_to_del = c(as.matrix(delete.top), as.matrix(delete.bottom))
    if (long == "TRUE") {
      return(list(curated, cells_to_del))
      
    } else {
      return(cells_to_del)
      
    }
  }
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")

#delete_cells(nk_plot[[2]],thr=50,bottom=4)
```
```{r make plots}
setwd("/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")

tumor_cells.plot  <- counts_plot(tumor_cells, ret = TRUE)
t_cells.plot      <- counts_plot(t_cells, ret = TRUE)
b_cells.plot      <- counts_plot(b_cells, ret = TRUE)
macrophages.plot  <- counts_plot(macrophages, ret = TRUE)
caf.plot          <- counts_plot(caf, ret = TRUE)
undef.plot        <- counts_plot(undef, ret = TRUE)
nk_plot           <- counts_plot(nk, ret = TRUE)

```
```{r deleting cells in all subsets}
#parameters chosen after graphs examination
require(dplyr)
nk.corr <- delete_cells(nk_plot[[2]],thr=50,top=4, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")

b_cells.cor <- delete_cells(b_cells.plot[[2]], thr=10, top=8 , bottom=1, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles", long=TRUE)

b_cells.cor.2 <- delete_cells(b_cells.cor[[1]], thr=28, top=4 , bottom=0, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")

caf.cor <- delete_cells(caf.plot[[2]], thr=10, top=15 , bottom=0, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")
caf.group <- caf.cor

t_cells.cor <- delete_cells(t_cells.plot[[2]],thr=15,top=10, bottom =0, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles", long=TRUE)

t_cells.cor.2 <- delete_cells(t_cells.cor[[1]], thr=3, top=0 , bottom=1, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")

macrophages.cor <- delete_cells(macrophages.plot[[2]],thr=12,top=2, bottom =3, path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles")
```

```{r create datasets for ICA}
#df - data frame of subset
#vec - cells to delete
#path - full path to destination
require(dplyr)

matrix_del_rows_scaled <- function (df, vec, path) {
  
  name <- deparse(substitute(df))
  
  df.rd <-
    data.frame(subset(df,!(rowname %in% vec))%>% select(-(tumor:cell_type)))
  
  rows <- df.rd[, 1]
  
  df.rd[, 1] <- NULL
  
  row.names(df.rd) <- rows
  
  df.scaled.t <- data.frame(t(scale(df.rd, scale = FALSE)))
  return(df.scaled.t)
}



nk.ica=matrix_del_rows_scaled(nk,nk.corr,path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/")
dim(nk.ica)
```

```{r doICABatch}
path_global="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/nk.ica_stability_2_25/"
name="nk.ica_stability_2_25"
vec=seq(2,25,1)

doICABatch <- function(df.scaled.t, vec, path_global) { 
  
  name <- deparse(substitute(df.scaled.t))
  print(name)
  ncomp=paste("[",paste(as.character(vec), collapse=" "),"]",sep=" ")
  
  setwd(path_global)
  dir.create(paste(name,"stability",vec[1],vec[length(vec)],sep="_"))
  print(paste(name,"stability",vec[1],vec[length(vec)],sep="_"))
  path_global=paste(path_global,paste(name,"stability",vec[1],vec[length(vec)],sep="_"),"/",sep="")
  print(path_global)
  setwd(path_global)
  
  write.table(
    df.scaled.t,
    file = paste(name, "_numerical.txt", sep = ""),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  fun="doICABatch"
  
  
  path=paste("'",path_global,"'",sep="")
  
  file=paste("'",paste(name, "_numerical.txt", sep = ""),"'",sep="")
  
  start="/Applications/MATLAB_R2013a_Student.app/bin/matlab -nodisplay -r \"cd('/Users/urszulaczerwinska/Documents/CURIE2015/ICA_Arnau/ICA_script/FastICA_25');"
  cmd=paste(start," ",fun,"(",path,",",file,",", ncomp,"); quit\"", sep="") 
  
  system(cmd)
  t.imp.path=paste(path_global,"avg_stability.plot.txt", sep="")
  T=read.delim(t.imp.path, header=FALSE, sep="\t")
  row.names(T)=as.character(T[,1])
  T=T[order(T[,1]),]
  png(filename=paste(path_global,"avg_stability.plot.png",sep=""))
  barplot(t(T[,2,drop=FALSE]), col='blue', ylim=c(0,1), main = paste(name,"stability",vec[1],vec[length(vec)],sep="_") )
  dev.off()
  
  return(T)
  # require(EBImage)
  # 
  # img = readImage(paste(path_global,"stability.png", sep="/"))
  # display(img, method = "raster")
}   


```

```{r doICA}

doICA<- function(df.scaled.t, path_global,n) {
  setwd(path_global)
  name <- deparse(substitute(df.scaled.t))
  dir.create(paste(name,n,sep="_"))
  
  path_global=paste(path_global,paste(name,n,sep="_"),"/",sep="")
  setwd(path_global)
  
  write.table(
    df.scaled.t,
    file = paste(name, "_numerical.txt", sep = ""),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  write.table(
    row.names(df.scaled.t),
    file = paste(name, "_ids.txt", sep = ""),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  write.table(
    colnames(df.scaled.t),
    file = paste(name, "_samples.txt", sep = ""),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  fun = "doICA"
  
  ncomp = n
  
  path = paste("'",path_global,"'",sep="")
  
  file = paste("'",paste(name, "_numerical.txt", sep = ""),"'",sep="")
  
  start = "/Applications/MATLAB_R2013a_Student.app/bin/matlab -nodisplay -r \"cd('/Users/urszulaczerwinska/Documents/CURIE2015/ICA_Arnau/ICA_script/FastICA_25');"
  cmd = paste(start," ",fun,"(",path,",",file,",", ncomp,"); quit\"", sep="") 
  
  system(cmd)
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
  
  return(list(A,S))
}   


```

```{r test for doICA}

path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/"
nk.doica=doICA(nk.ica, path , 3)      
nk.ica.batch=doICABatch(nk.ica, c(2,3,4), path)      

```
```{r testing combining scripts}
#df=nk.ica
running_batch_ica <- function (df, path) {
  name.batch <- paste(deparse(substitute(df)),"batch",sep=".")
  print(name.batch)
  #nk.ica.batch=doICABatch(nk.ica, seq(2,25,1), path)  
  assign(name.batch, doICABatch(df, seq(2,25,1), path))
  T <- get(name.batch)
  print(T)
  #n=T[which(T[,2]==max(T[,2])),1]
  #name.doica=paste(deparse(substitute(df)),"doica",sep=".")
  #print(name.doica)
  #assign(name.doica, doICA(df, path,n))
  #return(get(name.doica))
}
path="/Users/urszulaczerwinska/Documents/CURIE2015/Data/single.cell.melanoma/Average_Profiles/ICA_tirosh_subsets/"

running_batch_ica(nk.ica,path )

```

