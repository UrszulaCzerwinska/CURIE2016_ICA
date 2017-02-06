load("/bioinfo/users/uczerwin/SCDE_R/RDATA/estim_tirosh_light.RData")
library(scde)
estim_exprs <- function (df) {
  require(scde)
  df.t <- data.frame(t(df[, 5:ncol(df)]))
  colnames(df.t) <- df[, 1]
  df.t.unlog <-
    data.frame(apply(df.t, 2, function(x)
      as.integer(exp(x) - 1)))
  row.names(df.t.unlog) <- row.names(df.t)
  cd <- clean.counts(data.frame(as.matrix(df.t.unlog)))
  ifm <-
    scde.error.models(
      counts = cd,
      n.cores = 1,
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
                    n.cores = 1,
                    return.individual.posteriors = TRUE)
  return(list(cd=cd, ifm= ifm, prior=prior, posterior=posterior))
}

t1= Sys.time()
b_cells.estim <- estim_exprs(b_cells)
t2= Sys.time()
td=t1-t2 
save.image("/bioinfo/users/uczerwin/SCDE_R/RDATA/estim_tirosh_light110116.RData")


