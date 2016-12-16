path_global="/Users/ulala/Documents/CURIE/signatures/new_sig/"
setwd(path_global)
mcp_sets=data.frame(read.table(paste(path_global,"MCP",sep=""), header=TRUE, sep="\t"))
head(mcp_sets)
mcp_sets$id=1:nrow(mcp_sets)
mcp_sets$Cell.population = as.factor(mcp_sets$Cell.population)
w <- reshape(mcp_sets, 
             idvar = "id", 
             timevar = c("Cell.population"),
             drop="ENTREZID",
             direction = "wide")

# write.table(w,file="long_MCP.txt", sep="\t", quote=FALSE)
list_my=NULL
for (i in 2:ncol(w)) {
  
  list_i=na.omit(w[,i])[[1]]
  list_my=list(list_my,list_i)
  
  
}

l=apply(w, 2, function(x) print(array(t(na.omit(x)))))

paste(names(w), paste(l))
s=sapply(c(1:ncol(w)),function(i) paste(names(w)[i], paste(l[[i]],collapse = "    ") ,sep="    "))
write.table(s,file="long_MCP.txt", sep="\n", quote=FALSE, col.names=FALSE, row.names = FALSE) #pas terrible solution
      