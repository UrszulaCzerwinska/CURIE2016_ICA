
doICA<- function(df.scaled.t, path_global, n, name=FALSE, corr_folder="CORRELATION") {
  
  setwd(path_global)
  if (name != FALSE) {
    name <- name
  } else {
    name <- deparse(substitute(df.scaled.t))
  }
  
  dir.create(paste(name,n,sep="_"))
  
  path_global_1=paste(path_global,paste(name,n,sep="_"),"/",sep="")
  setwd(path_global_1)
  
  write.table(
    df.scaled.t,
    file = paste(name, "_numerical.txt", sep = ""),
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE,
    sep="\t"
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
  #ncomp = 2
  path = paste("'",path_global_1,"'",sep="")
  #path = "../../"
  file = paste("'",paste(name, "_numerical.txt", sep = ""),"'",sep="")
  #file ="myfile"
  #start = "/bioinfo/opt/build/Matlab-R2013a/bin/matlab -nodisplay -r \"cd('/data/users/uczerwin/ICA_Arnau/ICA_script/FastICA_25');"
  start = "/Applications/MATLAB_R2016a.app/Contents/MacOS/MATLAB_maci64 -nodisplay -r \"cd('/Users/ulala/Documents/CURIE/ICA_Arnau/ICA_script/FastICA_25');"
  time= paste("'",path_global_1,"/TimeSpent.txt'",sep="")
  #cmd = paste(start," ",fun,"(",path,",",file,",", ncomp,"); quit\"", sep="") 
  cmd = paste(start," tic; ",fun,"(",path,",",file,",", ncomp,"); TimeSpent = toc; dlmwrite(", time,",TimeSpent); quit\"", sep="") 
  
  system(cmd)
  system("rm *_numerical.txt")
}   

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
 
    path = "/Users/ulala/Documents/CURIE/Data/datasets_for_comparing_R_Matlab_fastica/UC/"
    setwd(path)
    df.doica = doICA(df.cen, path , 10, name = paste(dfr,"3",sep="_"))
    if (ncol(df) >= 50) {
      df.doica = doICA(df.cen, path , 50, name = paste(dfr,"3",sep="_"))
    }
    if (ncol(df) >= 100) {
      df.doica = doICA(df.cen, path , 100, name = paste(dfr,"3",sep="_"))
    }
  }
#}


run_ica_benchmark(dataset_list)