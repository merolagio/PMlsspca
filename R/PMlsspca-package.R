
#' Baseball hitters career and 1986 season total statistics
#' 
#' 16 statistics of 263 major league hitters some of the
#' overall career and others relative to the 1986 season. Available at StatLib.
#' The matrix has a block structure, defined by season offensive play, career 
#' offensive play and season defensive play.
#' 
#' 
#' @name hitters
#' @docType data
#' @format A \emph{263} by \emph{17} matrix. The first 16 variables are centered 
#' and standardized to unit norm.
#' \describe{
#' \item{TAB_86}{times at bat in 1986}
#' \item{HIT_86}{hits in 1986}
#' \item{HR_86}{home runs in 1986}
#' \item{RUN_86}{runs in 1986}
#' \item{RB_86}{runs batted-in in 1986}
#' \item{WAL_86}{walks in 1986}
#' \item{WAL}{walks during his career}
#' \item{PO_86}{put outs in 1986}
#' \item{ASS_86}{assists in 1986}
#' \item{ERR_86}{errors in 1986}
#' \item{YC}{years in the major leagues}
#' \item{TAB}{times at bat during his career}
#' \item{HIT}{hits during his career}
#' \item{HR}{home runs during his career}
#' \item{RUN}{runs during his career}
#' \item{RUNB}{runs batted-in during his career}
#' \item{Salary_86}{Salary in 1986}
#' }
#' @source \url{http://lib.stat.cmu.edu/datasets/baseball.data}
#' @keywords datasets
NULL


#' Baseball hitters statistics labels reference table
#' 
#' This data frame provides descriptive labels for the variables in the
#' hitters datasets matching the short ones used. 
#'   
#' @name hitters_labels
#' @docType data
#' @format A dataframe with columns
#' \describe{
#' \item{short.namethe}{labels ised in plots and tables}
#' \item{label}{the explanatory name}
#' \item{type}{whether offensive year 1986, defensive year 86 or offensive over the whole career}
#' }
#' @source \url{http://lib.stat.cmu.edu/datasets/baseball.data}
#' @keywords datasets
NULL


#' PMlspca: Least Squares Sparse PCA utilities
#'
#' PMlspca provides functions to compute Least Squares Sparse Principal Components Analysis
#' (LSSPCA) solutions, where sparsity is imposed while targeting PCA’s least-squares
#' reconstruction objective. This release accompanies the related article and is intended
#' to support full reproduction of the results reported therein.
#'
#' Computation relies on efficient C++ routines and includes multiple options for variable
#' selection and sparse loading estimation.
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
#' \code{\link{showload}} \tab Display `spca` loadings.\cr
#' \code{\link{compare.spca}} \tab Compare multiple `spca` objects (summaries and plots).\cr
#' \code{\link{aggregate_by_scale}} \tab Aggregate loadings/contributions by scale.\cr
#' }
#'
#' @docType package
#' @name PMlsspca-package
#' @import Rcpp, RcppEigen, ggplot2, RMTstats
#' @references Merola G. M. 2015 \emph{Sparse Principal Component Analysis: a Least Squares Approach}. Australian & New Zealand J. of Statistics. 57(3)\cr\cr
#' Merola G. M, Chen G. 2019 \emph{Projection sparse principal component analysis: An efficient least squares method}. Journal of Multivariate Analysis. 173. pp. 366-382.
#' @keywords package
#' @seealso \code{\link{lsspca}} for usage examples.
#' @useDynLib PMlsspca 
NULL
