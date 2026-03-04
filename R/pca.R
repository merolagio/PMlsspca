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
  #'######============================================================
  ## computes the PCA loadings with vexp, cumulative vexp and eigenvalues
######============================================================ 
  
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

#' \strong{Eigenvalue-based diagnostics for retaining principal components}
#'
#' Computes only the eigenvalues of a correlation or covariance matrix using
#' \code{EigenvaluesC()}, and provides standard diagnostics for selecting the
#' number of components to retain. The input \code{M} can be either a data matrix
#' (\eqn{n \times p}) or a symmetric correlation/covariance matrix (\eqn{p \times p}).
#'
#' Diagnostics include a screeplot, a Wachter (MP) QQ-plot, cumulative variance
#' explained (CVEXP), and (for correlation matrices) the Kaiser count (number of
#' eigenvalues larger than 1).
#'
#' @param M A numeric data matrix (\eqn{n \times p}) or a symmetric correlation /
#'   covariance matrix (\eqn{p \times p}). If \code{M} is square but not symmetric,
#'   it is treated as a data matrix and a warning is issued.
#' @param n Optional integer sample size. Required for the MP QQ-plot when
#'   \code{M} is already a correlation/covariance matrix. If \code{M} is a data
#'   matrix, \code{n} is set to \code{nrow(M)}.
#' @param make_cor Logical: only used if \code{M} is a data matrix. If
#'   \code{TRUE}, computes \code{cor(M)}; if \code{FALSE}, computes \code{cov(M)}.
#' @param nplot Integer: number of eigenvalues to include in the plots. Default
#'   is \code{p}.
#' @param kaiser_line Logical: if \code{TRUE}, adds a horizontal line at
#'   \eqn{y = 1} in the screeplot. Default is \code{TRUE} for correlation matrices
#'   and \code{FALSE} for covariance matrices.
#' @param nfit_line Integer: passed to \code{wachterqq()} to control the number
#'   of points used to fit the line.
#' @param rtn_scree Logical: should the screeplot object be returned?
#' @param prn_scree Logical: should the screeplot be plotted?
#' @param rtn_qq Logical: should the QQ-plot object be returned?
#' @param prn_qq Logical: should the QQ-plot be plotted?
#' @param rtn_values Logical: should eigenvalue statistics be returned?
#'   Default \code{TRUE}.
#' @param prn_cvexp integer: how many CVEXP values should be  printed as percentages (with a
#'   trailing \code{\%})?
#' @param digits_cvexp Integer: number of decimal digits used when printing
#'   CVEXP percentages.
#' @param tol_cor Numeric tolerance for detecting a correlation matrix from a
#'   symmetric square input \code{M} via \code{diag(M) == 1}.
#' @param tol_sym Numeric tolerance for checking symmetry of a square input
#'   matrix.
#'
#' @return A named list containing any combination of the following three
#'   elements, depending on the \code{rtn_*} flags. If no flag is \code{TRUE},
#'   returns \code{invisible(NULL)}.
#'
#' \describe{
#'   \item{\code{values}}{Returned if \code{rtn_values = TRUE}. A list with:
#'     \describe{
#'       \item{\code{eigvals}}{Eigenvalues in decreasing order.}
#'       \item{\code{vexp}}{Variance explained proportions (VEXP).}
#'       \item{\code{cvexp}}{Cumulative variance explained proportions (CVEXP).}
#'       \item{\code{kaiser}}{Number of eigenvalues larger than 1 for a
#'         correlation matrix; \code{NA} otherwise.}
#'       \item{\code{is_cor}}{Logical: whether \code{M} was treated as a
#'         correlation matrix.}
#'       \item{\code{n}}{Sample size used (or supplied).}
#'       \item{\code{p}}{Number of variables.}
#'     }
#'   }
#'   \item{\code{screeplot}}{Returned if \code{rtn_scree = TRUE}. The screeplot
#'     object produced by \code{screeplot()}; \code{NULL} otherwise.}
#'   \item{\code{wachterqq}}{Returned if \code{rtn_qq = TRUE}. The Wachter (MP)
#'     QQ-plot object produced by \code{wachterqq()}; \code{NULL} otherwise.}
#' }
#'
#' @details Eigenvalues are computed with the C++ routine \code{EigenvaluesC(M)}.
#' When \code{M} is a data matrix, \code{M} is replaced internally by
#' \code{cor(M)} or \code{cov(M)} depending on \code{make_cor}. Kaiser counts are
#' only meaningful for correlation matrices. The MP QQ-plot requires \code{n};
#' if \code{n} is missing when \code{M} is already a correlation/covariance
#' matrix, a warning is issued and the QQ-plot is skipped.
#'
#' @seealso \code{\link{screeplot}}, \code{\link{wachterqq}}, \code{\link{pca}}
#'
#' @export pc_retention
pc_retention <- function(M,
                         n            = NULL,
                         make_cor     = TRUE,
                         nplot        = NULL,
                         kaiser_line  = NULL,
                         nfit_line    = NULL,
                         rtn_scree    = FALSE, prn_scree    = TRUE,
                         rtn_qq       = FALSE, prn_qq       = TRUE,
                         rtn_values   = TRUE,  
                         prn_cvexp    = 0,
                         digits_cvexp = 1,
                         tol_cor      = 1e-8,
                         tol_sym      = 1e-8) {
  
  if (any(is.na(M)))
    stop("the matrix cannot contain NAs")
  if (!is.matrix(M)) M <- as.matrix(M)
  
  # --- Symmetry check first: decides data-matrix vs. square-matrix branch ---
  is_symmetric <- (nrow(M) == ncol(M)) &&
    (max(abs(M - t(M))) <= tol_sym)
  
  if (is_symmetric) {
    # M is already a cor/cov matrix - detect which one from the diagonal
    p      <- ncol(M)
    is_cor <- all(abs(diag(M) - 1) < tol_cor)
    
  } else {
    # M is a rectangular (or non-symmetric square) data matrix
    if (nrow(M) == ncol(M))
      warning("M is square but not symmetric - treating as a data matrix.")
    n  <- nrow(M)
    p  <- ncol(M)
    if (isTRUE(make_cor)) {
      M      <- stats::cor(M, use = "pairwise.complete.obs")
      is_cor <- TRUE
    } else {
      M      <- stats::cov(M, use = "pairwise.complete.obs")
      is_cor <- FALSE
    }
  }
  
  if (is.null(nplot))       nplot       <- p
  if (is.null(kaiser_line)) kaiser_line <- is_cor
  
  # --- Eigenvalues (C++) ---
  eigvals <- EigenvaluesC(M)
  eigvals <- sort(as.numeric(eigvals), decreasing = TRUE)
  
  # --- VEXP / CVEXP ---
  vexp  <- eigvals / sum(eigvals)
  cvexp <- cumsum(vexp)
  
  if (prn_cvexp > 0) {
    cvexp_chr        <- paste0(formatC(100 * cvexp[1:prn_cvexp], format = "f", digits = digits_cvexp), "%")
    names(cvexp_chr) <- paste0("PC", seq_along(cvexp_chr))
    print(cvexp_chr)
  }
  
  # Kaiser (only meaningful for correlation matrices)
  kaiser <- if (is_cor) sum(eigvals > 1) else NA_integer_
  
  # --- Plots ---
  scree_pl <- NULL
  qq_pl    <- NULL
  
  if (isTRUE(prn_scree) || isTRUE(rtn_scree)) {
    scree_pl <- screeplot(eigvals,
                          nplot       = nplot,
                          kaiser_line = isTRUE(kaiser_line),
                          ylab        = "eigenvalues",
                          prn         = isTRUE(prn_scree),
                          rtn         = isTRUE(rtn_scree))
  }
  
  if (isTRUE(prn_qq) || isTRUE(rtn_qq)) {
    if (is.null(n) || !is.finite(n) || n <= 1) {
      warning("To produce the MP QQ-plot you must supply `n` (sample size), unless `M` was a data matrix.")
      prn_qq <- FALSE
      rtn_qq <- FALSE
    } else {
      qq_pl <- wachterqq(eigvals   = eigvals,
                         p         = p,
                         n         = n,
                         cor       = is_cor,
                         nplot     = nplot,
                         nfit_line = nfit_line,
                         prn       = isTRUE(prn_qq),
                         rtn       = isTRUE(rtn_qq))
    }
  }
  
  # --- Assemble output: only include what was requested ---
  out <- list()
  
  if (isTRUE(rtn_values))
    out$values <- list(eigvals = eigvals,
                       vexp    = vexp,
                       cvexp   = cvexp,
                       kaiser  = kaiser,
                       is_cor  = is_cor,
                       n       = n,
                       p       = p)
  
  if (isTRUE(rtn_scree))
    out$screeplot <- scree_pl
  
  if (isTRUE(rtn_qq))
    out$wachterqq <- qq_pl
  
  if (length(out) == 0) invisible(NULL) else out
}
