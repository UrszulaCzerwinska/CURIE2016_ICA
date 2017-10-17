#####testing distribution for genes
#testing on NK

# nk.t <- data.frame(t(nk[,5:ncol(nk)]))
# dim(nk.t)
# head(nk.t)
# colnames(nk.t) <- nk[,1]
# nk.t.unlog <- data.frame(apply(nk.t, 2, function(x) as.integer(exp(x) - 1)))
# row.names(nk.t.unlog) <- row.names(nk.t)
# head(nk.t.unlog)

nk.t.unlog[1:10,1:10]
# LOAD LIBRARIES 
library(fitdistrplus) # fits distributions using maximum likelihood 
library(gamlss) # defines pdf, cdf of ZIP
# FIT DISTRIBUTION (mu = mean of poisson, sigma = P(X = 0) 
summary(t(nk.t.unlog[2,]))
i.vec =c(t(nk.t.unlog[17,]))
hist(i.vec)
quantile(i.vec)
plot(density(i.vec))

i.vec=c(0,63,1,4,1,44,2,2,1,0,1,0,0,0,0,1,0,0,3,0,0,2,0,0,0,0,0,2,0,0,0,0, 0,0,0,0,0,0,0,0,6,1,11,1,1,0,0,0,2)
mean(i.vec[i.vec > 0])
#fit_zip = fitdist(i.vec, 'ZIP', start = list(mu = 2, sigma = 0.5)) 
fit_zip = fitdist(i.vec, 'pois', start = list(lambda = 28))
#fit_zip = fitdist(i.vec, 'ZIP', start = list(mu = 2,  sigma = 0.5), lower=c(-Inf, 0.001), upper=c(Inf, 1), optim.method="L-BFGS-B")
fit_zip = fitdist(i.vec, 'gamma', start = list(shape =1, scale = 48), lower = c(0, 0))
#fit_zip = fitdist(i.vec, 'beta', start = list(shape1 =2, shape2 = 5), lower = c(0, 0), method = "mle")
fit_zip$estimate
fit_zip = fitdist(i.vec, 'nbinom', start = list(size=1, mu= 1), method = "mle")

rpois(1000, fit_zip$estimate)

# VISUALIZE TEST AND COMPUTE GOODNESS OF FIT 
plot(fit_zip) 
gofstat(fit_zip)

gene <- runif(rpois(1,fit_zip$estimate), min=0, max=125)
round(ppois(quantile(i.vec), fit_zip$estimate),2)
plot(density(rgamma(100, shape=fit_zip$estimate[1], scale= fit_zip$estimate[2])))
plot(density(rnbinom(100, size = fit_zip$estimate[1], mu = fit_zip$estimate[2])))


t_cells[1:10,1:10]
df=t_cells  
df.t <- data.frame(t(df[, 5:ncol(df)]))
colnames(df.t) <- df[, 1]
df.t.unlog <-
  data.frame(apply(df.t, 2, function(x)
    as.integer(exp(x) - 1)))
row.names(df.t.unlog) <- row.names(df.t)

i.vec =c(t(df.t.unlog[8,]))
hist(i.vec , prob=TRUE)            # prob=TRUE for probabilities not counts
lines(density(i.vec ))             # add a density estimate with defaults
lines(density(i.vec , adjust=2, kernel = "cosine", from= 0), lty="dotted", col= "blue", give.Rkern=TRUE)
d=density(i.vec , adjust=2, kernel = "cosine", from= 0)
hist(d$x*d$y)
dd <- approxfun(i.vec, method="constant")
hist(dd(1:10000))

hist(i.vec)
abline(v=mean(i.vec), lty=2)
points(mean(i.vec), dd(mean(i.vec)), cex=1.2, pch=20, col="blue")


hist(as.integer(d$x))
as.integer(i.vec)
mean(i.vec)
mean(as.integer(i.vec))
i.vec_plus= i.vec+0.00000001
fit_zip = fitdist(dd(1:10000), 'nbinom', start = list(size=0.3, mu= 24), lower = c(0, 0), control=list(trace=1, REPORT=1))
fit_zip = fitdist(i.vec, 'nbinom', start = list(size=0.3, mu= 2), lower = c(0, 0), control=list(trace=1, REPORT=1))

?dnbinom
data = rnbinom(1000, size = 1, mu = 4)
hist(data)
fit_zip = fitdist(as.integer(i.vec), 'nbinom', start = list(size=0.3, mu= 0.4), lower = c(0, 0), control=list(trace=1, REPORT=1))
fit_zip = fitdist(i.vec, 'nbinom', start = list(size=0.3, mu= 0.4), lower = c(0, 0), control=list(trace=1, REPORT=1))
fit_zip = fitdist(i.vec, 'ZIP', start = list(mu = 0.4, sigma = 1),  control=list(trace=1, REPORT=1), method= "mle") 
plot(fit_zip) 
gofstat(fit_zip)
fit_zip$estimate
write.table(t_cells.t, "/Users/ulala/Documents/CURIE/Data/single.cell.melanoma/Average_Profiles/t_cells_log.txt", sep="\t", quote = FALSE )


df.t.unlog=nk.t.unlog[1:100,]
i=17

####function 
estimate_negative_bin_parameters <- function(df) {
#df=caf
  require(fitdistrplus)
  require(gamlss)
  df.t <- data.frame(t(df[, 5:ncol(df)]))
  colnames(df.t) <- df[, 1]
  df.t.unlog <-
    data.frame(apply(df.t, 2, function(x)
      as.integer(exp(x) - 1)))
  row.names(df.t.unlog) <- row.names(df.t)
  
  my_params = rep(list(), nrow(df.t.unlog))
  for (i in 1:nrow(df.t.unlog)) {
    #i= 1967
    #print(i)
    i.vec = c(t(df.t.unlog[i,]))
    if ((length(which(i.vec > 0)) / length(i.vec) >= 0.25) || (is.na(i.vec)) ) {
      fit_zip = fitdist(i.vec,
                        'nbinom',
                        start = list(size = 0.3, mu = 24),
                        lower = c(0, 0))
      my_params[[i]] = fit_zip$estimate
    } else {
      my_params[[i]] = c(0, 0)
    }
    
  }
  return(my_params)
  
}
df.t.unlog[ 1967,]
estimate_negative_bin_parameters(macrophages)

#done:
#b_cels
#undef
#macrophages

#to do
#tumor


caf_nbin_params = estimate_negative_bin_parameters(caf)
t_cells_nbin_params = estimate_negative_bin_parameters(t_cells)
endothelial_cells_nbin_params = estimate_negative_bin_parameters(endothelial_cells)
tumor_cells_nbin_params = estimate_negative_bin_parameters(tumor_cells)
nk_nbin_params = estimate_negative_bin_parameters(nk)

#

j <- 10
sin_cell_data=caf_nbin_params
sampling_sing_cell <- function ( j, sin_cell_data) {
  i <- length(sin_cell_data)
  my_data <-  matrix(data = NA, nrow = i, ncol = j, byrow = TRUE)
  
  for (l in 1:i){
    my_data[l,] = rnbinom(j, size = sin_cell_data[[l]][[1]], mu = sin_cell_data[[l]][[2]])
    
    
  }
 # my_data[is.nan(my_data)] <- 0
  
  return(my_data)
}

hist(my_data[2,])

#simulate 500 caf cells 
cells_caf_sim <- sampling_sing_cell(500, caf_nbin_params )
df <- cells_caf_sim 

row.na <- apply(df, 1, function(X) any(is.nan(X)))
cells_caf_sim.non <- cells_caf_sim[!row.na,]
dim(cells_caf_sim.non)
#compute spectrum: eigen vectors and values

#mini function to transpose data
compute_trans_cell_type <- function(df) {
df.t <- data.frame(t(df[, 5:ncol(df)]))
colnames(df.t) <- df[, 1] 
return(df.t) }
#minifunction to unlog data
unlog_data <- function(df.t) {
  df.t.unlog <-
  data.frame(apply(df.t, 2, function(x)
    as.integer(exp(x) - 1)))
row.names(df.t.unlog) <- row.names(df.t)
return(df.t.unlog)
}
#covarance of real and simulated data without NaN rows
cov_cafs_sim <- cov(t(log1p(cells_caf_sim[!row.na,])))
cov_cafs <- cov(t(compute_trans_cell_type(caf)[!row.na,]))

dim(cov_cafs_sim )
dim(cov_cafs)

# compute_trans_cell_type(caf)[!row.na,][1:10,1:10]
# log1p(cells_caf_sim[!row.na,])[1:10,1:10]
# head(cov_cafs_sim[,1:10] )
# head(cov_cafs[,1:10])
#install package for big eigen values computation
#install.packages("RSpectra")
library("RSpectra")
#eigen of both
eigen_caf_sim <-eigs_sym(cov_cafs_sim, 100)
eigen_cafs <-eigs_sym(cov_cafs,100)

hist(eigen_cafs$values,100)
hist(eigen_caf_sim$values,100)

head(eigen_cafs$values)
head(eigen_caf_sim$values)

#try with tcells
#simulate 500 tcell cells 
cells_tcell_sim <- sampling_sing_cell(2040 , t_cells_nbin_params )
df <- cells_tcell_sim 

row.na <- apply(df, 1, function(X) any(is.nan(X)))

#cov
cov_tcell_sim <- cov(t(log1p(cells_tcell_sim [!row.na,])))
cov_tcell <- cov(t(compute_trans_cell_type(t_cells)[!row.na,]))

dim(cov_tcell_sim)
dim(cov_tcell)

#eigen of both
eigen_tcell_sim <-eigs_sym(cov_tcell_sim, 100)
eigen_tcell <-eigs_sym(cov_tcell,100)

barplot(eigen_tcell$values, main= "eigenvalues t_cells")
barplot(eigen_tcell_sim$values,  main= "eigenvalues t_cells simul")

eigen_tcell$vectors[,1]
eigen_tcell_sim$vectors[,1]
dist(rbind(eigen_tcell$vectors[,1],eigen_tcell_sim$vectors[,1]))^2
plot(eigen_tcell$vectors[,2], eigen_tcell_sim$vectors[,2], pch=16)
