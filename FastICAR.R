# install.packages("fastICA")
# fastICA(X, n.comp, alg.typ = c("parallel","deflation"),
#         fun = c("logcosh","exp"), alpha = 1.0, method = c("R","C"),
#         row.norm = FALSE, maxit = 200, tol = 1e-04, verbose = FALSE,
#         w.init = NULL)

library(fastICA)

# start.time <- Sys.time()
# 
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

time.res <- data.frame(system = character(), code = character(), ncomp = numeric(), df = character(), run = numeric(), time = numeric(),stringsAsFactors = FALSE)

counter = 0
for (ntimes in 1:3){
  
dataset_list = list("BRCATCGA", "METABRIC", "CHOL", "OVCA", "STAD_methylome")
list = dataset_list
#run_ica_benchmark <- funtion(list = list()) {
for (dfr in list) {

  df <-
    read.table(
      paste(
        "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/",
        dfr,
        ".txt",
        sep = ""
      ),
      sep = "\t",
      row.names = NULL,
      header = TRUE
    )
  row.names(df) <- df[, 1]
  df <- df[, 2:ncol(df)]
  #dim(df)
  if ("X" %in% colnames(df)) { df <- subset(df, select=-c(X)) }
  #dim(df)
  df.cen <- sweep(df, 1, rowMeans(df))
  df.cen.t <- t( df.cen)
  
  counter = counter +1
 
  
  start.time <- Sys.time()
  df.doica = fastICA(as.matrix(df.cen), 10, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 
  end.time <- Sys.time()
  time.taken <- end.time - start.time

  time.res[counter, ] <- c("R", "R", 10, dfr, ntimes, time.taken)
  
  if (ncol(df) >= 50) {
    counter = counter +1
    start.time <- Sys.time()
    df.doica = fastICA(as.matrix(df.cen), 50, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.res[counter, ] <- c("R", "R", 50, dfr, ntimes, time.taken)
    
  }
  if (ncol(df) >= 100) {
    counter = counter +1
   
    
    start.time <- Sys.time()
    df.doica = fastICA(as.matrix(df.cen), 100,,alg.typ = "deflation", maxit = 1000, verbose = TRUE) 
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    time.res[counter, ] <- c("R", "R", 100, dfr, ntimes, time.taken)
    
  }
  

}
}
time.res$time <- as.numeric(time.res$time)*60

write.table(time.res, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/R_time_benchmark_sec_deflation_tansposed.txt", sep="\t", quote=FALSE, row.names = FALSE)

start.time <- Sys.time()

df.doica = fastICA(as.matrix(df.cen), 100,alg.typ = "deflation", maxit = 1000, verbose = TRUE) 
end.time <- Sys.time()
time.taken <- end.time - start.time



time_all <- read.table(
  
  "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/time_benchmark_all.txt",
  sep = "\t",
  row.names = NULL,
  header = TRUE
)

time_defl <-subset(time_all, type == "deflation")

(p1 <- ggplot(data = time_defl , aes(x = ncomp, y = time, colour = df)) +       
    geom_line() + geom_point(aes(shape=code)))

time_R <-subset(time_defl, code == "MATLAB")

(p2 <- ggplot(data = time_R , aes(x = ncomp, y = time, colour = df)) +       
    geom_line() + geom_point())

time_paral <-subset(time_all, type == "parallel")

time_MARLAB <-subset(time_defl, code == "MATLAB")

(p3 <- ggplot(data = time_paral , aes(x = ncomp, y = time, colour = df)) +       
    geom_line() + geom_point())

library(ggplot2)


#redo and save METABRIC to compare

dfr <- "METABRIC"
df <-
  read.table(
    paste(
      "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/",
      dfr,
      ".txt",
      sep = ""
    ),
    sep = "\t",
    row.names = NULL,
    header = TRUE
  )
row.names(df) <- df[, 1]
df <- df[, 2:ncol(df)]


#dim(df)
if ("X" %in% colnames(df)) { df <- subset(df, select=-c(X)) }
#dim(df)
df.cen <- sweep(df, 1, rowMeans(df))
METABRIC.doica.def_10 = fastICA(as.matrix(df.cen), 10, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 
METABRIC.doica.def_50 = fastICA(as.matrix(df.cen), 50, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 
METABRIC.doica.def_100 = fastICA(as.matrix(df.cen), 100, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 

write.table(METABRIC.doica.def_10$S, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/S_METABRIC_R_def_10.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(t(METABRIC.doica.def_10$A), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/A_METABRIC_R_def_10.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(METABRIC.doica.def_50$S, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/S_METABRIC_R_def_50.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(t(METABRIC.doica.def_50$A), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/A_METABRIC_R_def_50.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(METABRIC.doica.def_100$S, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/S_METABRIC_R_def_100.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(t(METABRIC.doica.def_100$A), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/A_METABRIC_R_def_100.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

write.table(data.frame( GENES =row.names(df) , METABRIC.doica.def_10$S), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/S_METABRIC_R_def_10.num_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.frame( GENES =row.names(df) ,METABRIC.doica.def_50$S), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/S_METABRIC_R_def_50.num_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.frame( GENES =row.names(df) ,METABRIC.doica.def_100$S), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/S_METABRIC_R_def_100.num_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)




METABRIC.doica.paral_10 = fastICA(as.matrix(df.cen), 10, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 
METABRIC.doica.paral_50 = fastICA(as.matrix(df.cen), 50, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 
METABRIC.doica.paral_100 = fastICA(as.matrix(df.cen), 100, alg.typ = "deflation", maxit = 1000,  verbose = TRUE) 

write.table(METABRIC.doica.paral_10$S, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/S_METABRIC_R_paral_10.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(t(METABRIC.doica.paral_10$A), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/A_METABRIC_R_paral_10.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(METABRIC.doica.paral_50$S, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/S_METABRIC_R_paral_50.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(t(METABRIC.doica.paral_50$A), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/A_METABRIC_R_paral_50.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(METABRIC.doica.paral_100$S, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/S_METABRIC_R_paral_100.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(t(METABRIC.doica.paral_100$A), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/A_METABRIC_R_paral_100.num", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

write.table(data.frame( GENES =row.names(df) , METABRIC.doica.paral_10$S), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/S_METABRIC_R_paral_10.num_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.frame( GENES =row.names(df) ,METABRIC.doica.paral_50$S), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/S_METABRIC_R_paral_50.num_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

write.table(data.frame( GENES =row.names(df) ,METABRIC.doica.paral_100$S), file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/S_METABRIC_R_paral_100.num_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



##########
path_file <- "/Users/ulala/Documents/CURIE/Data/Tatiana_TCGA/M7_CELLCYCLE.rnk"
path_file <- "/Users/ulala/Documents/CURIE/Data/Tatiana_TCGA/M3_SMOOTH_MUSCLE.rnk"
path_file <- "/Users/ulala/Documents/CURIE/Data/Tatiana_TCGA/M12_MYOFIBROBLASTS.rnk"

path_file <- "/Users/ulala/Documents/CURIE/Data/Tatiana_TCGA/immune_metagene.txt"




  
  IC = read.table(path_file, sep= "\t", stringsAsFactors = FALSE, header = FALSE)
  CELL_CYCLE = IC[order(IC[,1]),]
  SMOOTH_MUSCLE = IC[order(IC[,1]),]
  MYOFIBROBLASTS = IC[order(IC[,1]),]
  IMMUNE = IC[order(IC[,1]),]
  ADIPOSE = IC[order(IC[,1]),]
  common_genes <- intersect(intersect(intersect(CELL_CYCLE[,1],SMOOTH_MUSCLE[,1]  ),MYOFIBROBLASTS[,1]), IMMUNE[,1])
  
  Biton_metagenes <- data.frame(GENES = CELL_CYCLE[which( CELL_CYCLE[,1] %in% common_genes ),1], CELL_CYCLE = CELL_CYCLE[which( CELL_CYCLE[,1] %in% common_genes ),2],  SMOOTH_MUSCLE = SMOOTH_MUSCLE[which(SMOOTH_MUSCLE[,1] %in% common_genes ),2],MYOFIBROBLASTS = MYOFIBROBLASTS[which(MYOFIBROBLASTS[,1] %in% common_genes ),2], IMMUNE = IMMUNE[which( IMMUNE[,1] %in% common_genes ),2] )
  
  
  write.table( Biton_metagenes, file= "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/CORRELATION/Biton_for_corel_S.xls", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  
