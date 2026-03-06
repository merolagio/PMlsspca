#' \strong{Computes principal components solutions}
#' 
#' Creates PCA output as an \link{PMlsspca-package} object
#' which can run on the spca utilities. .
#' 
#' @param M A data matrix or correlation or covariance matrix.
#' @param ncomps Integer: number of loadings to retain. If missing all loadings are retained.
#' @param centerdata = FALSE, should variables be centered to zero mean automatically done if any mean is not zero
#' @param scaledata = FALSE , should variables be centered to unit variance?
#' @param screeplot Logical: should the screeplot be plotted?
#' @param kaiser_line Logical: adds a horizontal line to the screeplot at y = 1 
#' @param kaiser.print Logical: should the kaiser rule be computed, printed and returned?.
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
#' \item{eigenvalues}{all the eigenvalues}
#' \item{kaiser}{The number of eigenvalues larger than one, if \code{kaiser.print = TRUE.}}
#' \item{method}{Always PCA}
#' \item{Call}{The called with its arguments.}
#' }
#' @details \emph{ncomps} is just the number of components retained, doesn't speed up the function. Kaiser rule determines the number of components as the number of eigenvalues larger than one. It should be used only for correlation matrices, if called on a covariance matrix a warning is generated. 
# 
#' @export pca
#' @seealso See also \code{\link{print.spca}, \link{summary.spca}}
pca <- function(M, ncomps, centerdata = FALSE, scaledata = FALSE,
                 screeplot= FALSE, kaiser_line = F, kaiser.print = FALSE){
  ######============================================================P
  ## computes the PCA loadings with vexp, cumulative vexp and eigenvalues
  ######============================================================P 
  
  if(any(is.na(M)))
    stop("The matrix cannot contain NAs")
  
  if (is.data.frame(M))
    M = as.matrix(M)
  if (!is.matrix(M))
    stop("M must be a matrix")
  
  n = nrow(M)
  p = ncol(M)
  
  
  if(missing(ncomps))
    ncomps = p
  if((ncomps > p) || (ncomps <= 0)){
    warning(paste("incorrect value for ncomps in pca, set to), p"))
    ncomps = p}
  
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

  
    ee$vec = scaleColsC(ee$vec[, 1:ncomps], 0, sign(ee$vec[1, 1:ncomps]))
    rownames(ee$vec) = colnames(S)
    colnames(ee$vec) = paste0("PC", 1:ncomps)
    contributions = scaleColsC(ee$vec[,1:ncomps], 1, rep(1, ncomps))
    rownames(contributions) = colnames(S)
    ldlst = lapply(1:ncomps, function(i, x) x[, i], x = ee$vec[, 1:ncomps])
    vexp = ee$val[1:ncomps]/sum(ee$val)
  
  out = list(loadings = ee$vec[, 1:ncomps, drop = FALSE], 
             contributions = contributions,
             ncomps = ncomps, cardinality = rep(p, ncomps), 
             ind = as.list(rep(list(1:p), ncomps)), 
             vexp = vexp, vexpPC = vexp, cvexp = cumsum(ee$val[1:ncomps])/sum(ee$val),
             loadlist = ldlst, eigenvalues = ee$val)
  if (is_datamatrix_M){
    out$scores = M %*% out$loadings[, 1:ncomps, drop = FALSE]
    colnames(out$scores) = paste0("Comp", 1:ncomps)
    
  if (out$ncomps > 1){
      out$corComp = MkCorCompMat(out$loadings, S, d = out$ncomps)
  }
  }
  class(out) = c("spca", "list")
  out$vif <- make_vif_R(out, S, prn = FALSE)
  
  if (screeplot == TRUE){
    pl = screeplot(ee$val, kaiser_line = kaiser_line)
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
}
