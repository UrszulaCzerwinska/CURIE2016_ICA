unit_vec <- function(x) {
  unit_vecs <- apply(x, 2, function(vec) vec / sqrt(sum(vec ^ 2)))
  return(unit_vecs)
}

S_unit <- unit_vec(S_TcellLiver_numerical.txt_12)
dim(S_unit)
dim(S_TcellLiver_numerical.txt_12)
res <- t(S_unit)%*%as.matrix(TcellLiver_numerical)
colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}
rowVars<-function (x,na.rm = TRUE)
{
  sqr = function(x) x * x
  n = rowSums(!is.na(x))
  n[n <= 1] = NA
  return(rowSums(sqr(x - rowMeans(x,na.rm = na.rm)), na.rm = na.rm)/(n - 1))
}

hist(rowVars(res))
A_unit <- unit_vec(A_TcellLiver_numerical.txt_12)
res2 <- t(A_unit)%*%as.matrix(t(TcellLiver_numerical))
hist(rowVars(res2))
length(rowVars(res2))
fev1 <- rowVars(res2)/sum(var(t(TcellLiver_numerical)))
fev2 <- rowVars(res)/sum(var(t(TcellLiver_numerical)))
hist(fev1*nrow(TcellLiver_numerical))
hist(fev2*nrow(TcellLiver_numerical))
