require(scde)
dim(nk)
nk.t <- data.frame(t(nk[,5:ncol(nk)]))
dim(nk.t)
head(nk.t)
colnames(nk.t) <- nk[,1]
nk.t.unlog <- data.frame(apply(nk.t, 2, function(x) as.integer(exp(x) - 1)))
row.names(nk.t.unlog) <- row.names(nk.t)
head(nk.t.unlog)
cd.nk <- clean.counts(data.frame(as.matrix(nk.t.unlog)))

# > dim(cd.nk)
# [1] 7286   50
# > dim(nk.t)
# [1] 23686    51

o.ifm.nk <- scde.error.models(counts = cd.nk, n.cores = 4, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1, min.size.entries	 =5)
valid.cells.nk <- o.ifm.nk$corr.a > 0
table(valid.cells.nk)
o.ifm.nk <- o.ifm.nk[valid.cells.nk, ]
dim(o.ifm.nk)
prior.nk <- scde.expression.prior(models = o.ifm.nk, counts = cd.nk, length.out = 400, show.plot = FALSE)
head(prior.nk)
posterior.nk <- scde.posteriors(o.ifm.nk, cd.nk, prior.nk, n.cores = 4,return.individual.posteriors =TRUE)


##########################
df=undef

estim_exprs <- function (df, min.lib.size=2300) {
  require(scde)
  df.t <- data.frame(t(df[, 5:ncol(df)]))
  colnames(df.t) <- df[, 1]
  df.t.unlog <-
    data.frame(apply(df.t, 2, function(x)
      as.integer(exp(x) - 1)))
  row.names(df.t.unlog) <- row.names(df.t)
  cd <- clean.counts(data.frame(as.matrix(df.t.unlog)), min.lib.size = min.lib.size)
  ifm <-
    scde.error.models(
      counts = cd,
      n.cores = 4,
      threshold.segmentation = TRUE,
      save.crossfit.plots = FALSE,
      save.model.plots = FALSE,
      verbose = 1,
      min.size.entries	 = 5
    )
  valid.cells <- ifm$corr.a > 0
  ifm <- ifm[valid.cells,]
  prior <-
    scde.expression.prior(
      models = ifm,
      counts = cd,
      length.out = 400,
      show.plot = FALSE
    )
  posterior <-
    scde.posteriors(ifm,
                    cd,
                    prior,
                    n.cores = 4,
                    return.individual.posteriors = TRUE)
  return(list(cd=cd, ifm= ifm, prior=prior, posterior=posterior))
}



########################
plot_posterior(posterior.nk, 3, prior.nk)
title(paste(rownames(cd.nk)[3]))
express <- prior.nk$x
prob <- posterior.nk$jp[3, ]

p.nk = 0.2
N = 10000
val.3 <- sample(express, size = p.nk * N, replace = TRUE, prob = prob)
log(sum(exp(val.3)))
plot(density(val.3))

nk.exprs <- data.frame(rep(0, nrow(nk.t)))
row.names(nk.exprs) <- row.names(nk.t)


#### quite fast - computes sum of expression for a given nr of single cells for nk
t1= Sys.time()
for (i in 1:nrow(cd.nk)) {
  row.name.nk <- row.names(cd.nk)[i]

  express <- prior.nk$x
  prob <- posterior.nk$jp[i,]
  val.n <-
    sample(express,
           size = p.nk * N,
           replace = TRUE,
           prob = prob)
  nk.exprs[row.name.nk,] <- sum(exp(val.n))
}

length(which(nk.exprs > 0))
t2= Sys.time()
t1-t2 #Time difference of -5.87026 secs

#########################################


caf.estim <- estim_exprs(caf)
macro.estim <- estim_exprs(macrophages)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh060117#2.RData")
t_cell.estim <- estim_exprs(t_cells)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh130117#4.RData")
b_cell.estim <- estim_exprs(b_cells)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh160117#1.RData")

str(t_cell.estim)
dim(cd)

#missing endothelial & tumor & undef
b_cell.estim <- estim_exprs(b_cells)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh160117#1.RData")
endothelial_cells.estim <- estim_exprs(endothelial_cells)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh230117#1.RData")
tumor_cells.estim <- estim_exprs(tumor_cells, 7000)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh_tumor.RData")
undef.estim <- estim_exprs(undef,2300)
save.image("~/Documents/CURIE/R_code/RDATA/estim_tirosh240117#2.RData")
