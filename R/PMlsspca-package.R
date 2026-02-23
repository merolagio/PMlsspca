
#' PMlsspca: Least Squares Sparse PCA utilities
#'
#' Tools to compute LSSPCA solutions and reproduce the results reported in the accompanying article.
#'
#' The package provides functions to compute Least Squares Sparse Principal Components Analysis (LSSPCA) solutions, where sparsity is imposed while targeting PCA’s least-squares reconstruction objective. This release accompanies the related article and is intended  to support full reproduction of the results reported therein.
#'
#' Computation relies on efficient C++ routines and includes multiple options for variable selection and sparse loading estimation.
#'
#' S3 methods for objects of class `spca` include:
#' \tabular{ll}{
#' \code{\link{print.spca}} \tab Print nonzero loadings.\cr
#' \code{\link{plot.spca}} \tab Plot variance explained and sparse loadings.\cr
#' \code{\link{summary.spca}} \tab Summarize the solution with key statistics.\cr
#' }
#'
#' Additional utilities include:
#' \tabular{ll}{
#' \code{\link{showload}} \tab prints, the nonzero loadings or contributions separatelyfor each sPC.\cr
#' \code{\link{compare.spca}} \tab Compare multiple `spca` objects (summaries and plots).\cr
#' \code{\link{aggregate_by_scale}} \tab Aggregate loadings/contributions by scale.\cr
#' }
#'
#' The script to replicate the examples in the article is available at
#' \code{inst/scripts/Replicate_Examples.R} and can be located after installation via \code{system.file("scripts", "Replicate_Examples.R", package = "PMlsspca")}.
#' @keywords internal
"_PACKAGE"
#' @useDynLib PMlsspca, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Sparse principal components object
#'
#' Objects of class \code{"spca"} are returned by \code{\link{lsspca}} (and related functions). They store sparse loadings, scores, and summary quantities for the selected components.
#'
#' @section Components:
#' An object of class \code{spca} is a list with the following elements:
#'  
#' \describe{
#' Essential
#' \item{\emph{loadings}}{Numeric matrix \eqn{p \times r} of sparse loadings.}
#' \item{\emph{vexp}}{Vector of variance explained by each sparse component.}
#' \item{\emph{vexpPC}}{a vector of variance explained by the PCss named.} 
#' \item{\emph{ncomps}}{an integer declaring the number of components named.}
#' Optional
#' \item{\code{contributions}}{Numeric matrix \eqn{p \times r} of contributions (sign convention as stored).}
#' \item{\code{cvvexp}}{a vector of cumulative variance explained.}
#' \item{\code{rpcvexp}}{a vector of recovered variance relative to the corresponding PC.}
#' \item{\code{scores}}{Numeric matrix of component scores. Returned only if the data matrix is passed.}
#' \item{\code{card}}{Optional integer vector of component cardinalities (number of nonzeros).}
#' \item{\code{loadingslist}}{Optional list of per-component sparse loading vectors.}
#' \item{\code{corComp}}{Optional \eqn{r \times r} correlation matrix among sparse components.}
#' }
#'
#' @seealso \code{\link{print.spca}}, \code{\link{plot.spca}}, \code{\link{summary.spca}}
#'
#' @name spca_object
NULL
