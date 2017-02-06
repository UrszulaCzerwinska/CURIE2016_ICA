---
title: "Tirosh et al. distribution estimation"
output: html_notebook
---

```{r, message=FALSE, warning=FALSE, include=FALSE}
install.packages("pscl")
install.packages("MASS")
install.packages("devtools")

require(devtools)
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)

```

```{r}
tumor_cells[1:10,1:10]
dim(tumor_cells)
tc.num.test1 <- data.frame(tumor_cells[1:10,5:50])
head(tc.num.test1)
tc.t.test1 <- data.frame(t(tc.num.test1))
data(es.mef.small)
head(es.mef.small)
str(es.mef.small)

colnames(tc.t.test1 ) <- tumor_cells[1:10,1]
rownames(tc.t.test1 ) <- colnames(tumor_cells)[5:50]
cd <- clean.counts(tc.num.test1, min.lib.size=500, min.reads = 1, min.detected = 1)
cd1 <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
str(es.mef.small)
summary(es.mef.small)
dim(es.mef.small)
str(tc.num.test1)
summary(tc.num.test1[,1:10])
dim(tc.num.test1)
#tc.t.test1 <- data.frame(t(tc.num.test1))
summary(tc.t.test1)

tc.num.test1.unlog=data.frame(apply(tc.t.test1,2,function(x) as.integer(exp(x) - 1)))
#as.integer(exp(tc.num.test1[,2]) -1)

summary(tc.num.test1.unlog)
str(tc.num.test1.unlog)

cd2 <- clean.counts(data.frame(as.matrix(tc.num.test1.unlog)), min.lib.size=1, min.reads = 1, min.detected = 1)
cd3 <- clean.counts(data.frame(as.matrix(tc.t.test1 )), min.lib.size=1, min.reads = 1, min.detected = 1)

o.ifm <- scde.error.models(counts = cd3, n.cores = 4, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1,min.size.entries	 =5)
head(o.ifm)
rownames(o.ifm)
o.prior <- scde.expression.prior(models = o.ifm, counts = cd2, length.out = 400, show.plot = TRUE)
ediff <- scde.expression.difference(o.ifm, cd2, groups=groups, o.prior, n.randomizations  =  100, n.cores  =  4, verbose  =  1, return.posteriors = TRUE)
groups <- as.factor(c(rep("1", 1), rep("2",9)))
ediff$joint.posteriors[[1]]

pst <- scde.posteriors(o.ifm, cd2, o.prior, n.cores = 4)
max(o.prior$x)
plot(y = pst[9,], x = seq(0,6, length.out=401), type= "l")
```

```{r}

require(scde)
data(es.mef.small)

# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
# the group factor should be named accordingly
names(sg) <- colnames(es.mef.small)  
table(sg)
cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 4, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
data(o.ifm)
head(o.ifm)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
# define two groups of cells
groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels  =  c("ESC", "MEF"))
names(groups) <- row.names(o.ifm)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)
scde.browse.diffexp(ediff, o.ifm, cd, o.prior, groups = groups, name = "diffexp1", port = 1299)
scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior)

```


