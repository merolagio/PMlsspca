#' \strong{Computes principal components solutions}
#' 
#' Useful for creating an object with stanrd PCA loadings of class \code{spca}
#' which can run on the spca utilities. Other more comprehensive functions
#' for PCA exist in R, such as \code{princomp}, for example.
#' 
#' @param M A data matrix or correlation or covariance matrix.
#' @param ncomps Integer: number of loadings to retain. If missing all loadings are
#' retained.
#' @param only.values Logical: should only the eigenvalues be computed?
#' @param screeplot Logical: should the screeplot be plotted?
#' @param kaiser.print Logical: should the kaiser rule be computed, printed and
#' returned?.
#' @return An object of class \emph{spca} is returned, which contains:
#' \describe{
#' \item{loadings}{Matrix with the loadings scaled to unit \eqn{L_2} norm in 
#' the columns.} 
#' \item{contributions}{Matrix of loadings scaled to unit \eqn{L_1} norm.} 
#' \item{ncomps}{integer number of components computed.} 
#' \item{cardinality}{Vector with the cardinalities of each loadings.} 
#' \item{ind}{List with the indices of the non-zero loadings for each component.} 
#' \item{vexp}{Vector with the \% variance explained by each component.} 
#' \item{cvexp}{Vector with the \% cumulative variance explained by each component.} 
#' \item{vexpPC}{Vector with the \% variance explained by each principal component.} 
#' \item{rcvexp}{Vector with the proportion of cumulative variance explained by each 
#' component over the cumulative variance explained by the corresponding PCs.}
#' \item{Call}{The called with its arguments.}
#' \item{kaiser}{The number of eigenvalues larger than one, if \code{kaiser.print = TRUE}.}
#' }
#' @details \emph{ncomps} is just the number of components retained from the full eigen
#' decomposition, doesn't speed up the function. \emph{only.values} does not
#' compute the loadings and is more efficient. Kaiser rule determines the
#' number of components as the number of eigenvalues larger than one. It should
#' be used only for correlation matrices, if called on a covariance matrix a
#' warning is generated.
#' @examples
#'  \dontrun{ 
#' 
#' }
#' @export pca
#' @seealso See also \code{\link{print.spca}, \link{summary.spca}}
pca <- function(M, ncomps, centerdata = FALSE, scaledata = FALSE,
                only.values = FALSE, screeplot= FALSE, 
                kaiser.print = FALSE){
  #'######=============================================================  
  ## computes the PCA loadings with vexp, cumulative vexp and eigenvalues
  #'######=============================================================  
  
  if (is.data.frame(M))
    M = as.matrix(M)
  if (!is.matrix(M))
    stop("M must be a matrix")
  
  n = nrow(M)
  p = ncol(M)
  is_datamatrix_M = FALSE
  if ((n != p) & !isSymmetric(M))
    is_datamatrix_M = TRUE
  if (is_datamatrix_M == TRUE){
   message("input is data matrix, computing covariance matrix")
    if (centerdata | scaledata){
      M = scaleC(M, centerdata, scaledata)
     }
      S = crossprod(M)
  } 
  else{
    S = M
  }

  ee = EigenC(as.matrix(S))

  if(missing(ncomps))
    ncomps = p
  if (only.values == FALSE){
    ee$vec = scaleColsC(ee$vec[, 1:ncomps], 0, sign(ee$vec[1, 1:ncomps]))
    rownames(ee$vec) = colnames(S)
    colnames(ee$vec) = paste0("PC", 1:ncomps)
  }
  else
    ee$vec = NULL
  vexp = ee$val[1:ncomps]/sum(ee$val)
  contributions = scaleColsC(ee$vec[,1:ncomps], 1, rep(1, ncomps))
  rownames(contributions) = colnames(S)
  out = list(loadings = ee$vec[, 1:ncomps, drop = FALSE], 
             contributions = contributions,
             ncomps = ncomps, cardinality = rep(p, ncomps), 
             ind = as.list(rep(list(1:p), ncomps)), 
             vexp = vexp, vexpPC = vexp, cvexp = cumsum(ee$val[1:ncomps])/sum(ee$val),
             loadlist = lapply(1:ncomps, function(i, x) x[, i], x = ee$vec))
  
  if (is_datamatrix_M == TRUE){
    out$scores = M %*% out$loadings[, 1:ncomps, drop = FALSE]
    colnames(out$scores) = paste0("Comp", 1:ncomps)
  }  
  if (out$ncomps > 1){
      out$corComp = MkCorCompMat(out$loadings, S, d = out$ncomps)
  }
  class(out) = c("spca", "list")
  out$vif <- make_vif_R(out, S, prn = FALSE)
  
  old.par <- par(no.readonly = TRUE) 

  if (screeplot == TRUE){
    toplo = 100*ee$val/sum(ee$val)
    par(mar = c(4, 4, 1, 1))
    plot(1:p, toplo, xlab = "component", ylab = "% eigenvalue", type = "b", pch = 16, xaxt = "n")
    axis(side = 1, at = 1:p)
    if (kaiser.print)
      abline(a = 1, b = 0, lty = 2)
  }
  if (kaiser.print){
    if (any((diag(S) - 1) > 1e-6))
      warning("kaiser rule should be used only for correlation matrices")
    kaiser =   sum(ee$val > 1.0)
    out$kaiser = kaiser
    print(paste("number of eigenvalues larger than 1 is", kaiser ))
  }
  out$method = "PCA"
  out$Call = match.call()
  return(out)
  on.exit(par(old.par))
}
