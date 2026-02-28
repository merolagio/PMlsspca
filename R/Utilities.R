

#' Replication script for the article
#'
#' Returns the path to the script that reproduces the examples reported in the article.
#'
#' @return A character string giving the path to `Replicate_Examples.R` in the installed package (inst/scripts/"Replicate_Examples.R).
#' @export
replicate_article_examples <- function() {
  system.file("scripts", "Replicate_Examples.R", package = "PMlsspca")
}

list2vec = function(li, uniq = T, sorted = F){
  a = c(unlist(li))
  if(uniq) a = unique(a)
  if(sorted)
    a = sort(a)
  return(a)
}

vec2list = function(vec){
  u = unique(vec)
  li = vector("list", length(u))
  for (j in 1:length(u))
    li[[j]] = which(vec == u[j])
  names(li) = u
  li
} 

fac2list = function(fac){
  u = levels(fac)
  li = vector("list", length(u))
  for (j in 1:length(u))
    li[[j]] = which(fac == u[j])
  names(li) = u
  li
} 

list2fac = function(x){
  ca = sapply(x, length)  
  if(is.null(names(x)) )
    namex = paste0("e", 1:length(x))
  else namex = names(x)
  fa = c(unlist(x))
  fa = rep(0, length(fa))
  for (i in 1:length(x)){
    fa[x[[i]]] = i
  }
  return(factor(fa, labels = namex))
}

vec2fac = function(v){
  u = unique(v)
  u
  val = rep(0, length(v))
  for(i in 1:length(u)){
    val[v == u[i]] = i
  }
  factor(val, labels = u)
}
#' Wachter (Marchenko--Pastur) QQ plot for eigenvalues
#'
#' Produces a QQ-plot comparing observed eigenvalues to Marchenko--Pastur
#' (Wachter) theoretical quantiles for aspect ratio \eqn{\gamma = n/p}.
#'
#' @param eigvals Numeric vector of eigenvalues (assumed sorted decreasing).
#' @param p Integer. Number of variables.
#' @param n Integer. Sample size.
#' @param gamma Numeric. Aspect ratio; defaults to `n/p` if missing.
#' @param cor Logical. If `TRUE`, rescales MP quantiles to match a correlation
#'   matrix trace (sum to `p`).
#' @param nplot Integer. Number of leading eigenvalues to include; defaults to
#'   `length(eigvals)`.
#' @param nfit_line Integer or `NULL`. If positive, fits an `lm` line using the
#'   last `nfit_line` points; if negative, excludes the nfit largest values.
#' @param addtitle Logical. If `TRUE`, adds a plot title.
#' @param prn Logical. If `TRUE`, prints the plot.
#' @param rtn Logical. If `TRUE`, returns the ggplot object.
#'
#' @return If `rtn = TRUE`, a `ggplot` object; otherwise `NULL` (invisibly).
#'
#' @examples
#' # wachterqq(eigvals, p = ncol(X), n = nrow(X), cor = TRUE, nfit_line = 5, rtn = TRUE)
#'
#' @export
wachterqq = function(eigvals, p, n, gamma, cor = T, nplot, nfit_line = NULL, addtitle = TRUE, prn = TRUE, rtn = FALSE){
  
  if(missing(gamma)) gamma = n/p
  
  if(missing(nplot)) nplot = length(eigvals)
  
  probs <- ((p - (1:p) + 1)- 0.5)/p
  mp_quantiles <- RMTstat::qmp(p = probs, svr = gamma)
  
  if(cor) mp_quantiles <- p * mp_quantiles/sum(mp_quantiles)
  
  df = data.frame(expected = mp_quantiles[1:nplot], observed = eigvals[1:nplot])
  pl = ggplot(df, aes(x = expected, y = observed)) + geom_point(size = 2) + theme_tuto()
  
  if ((is.numeric(nfit_line)) && (nfit_line != 0)){
    if (nfit_line < 0) nfit_line = nplot + nfit_line 
    pl = pl + geom_smooth(data = df[(nplot - nfit_line + 1):nplot, ], se = F, method = "lm")
    }
  if(addtitle)
    pl = pl + labs(title = "wachter qq-plot") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  if(prn)
    print(pl)
  if(rtn)
    return(pl)
}

#' \strong{Scree plot of eigenvalues}
#'
#' Plots the first `nplot` eigenvalues (or their proportions) against component order.
#'
#' @param x A numeric vector of eigenvalues, or a list containing a numeric element named `values`.
#' @param nplot Integer. Number of leading eigenvalues to plot; defaults to `length(x)` (or `length(x$values)` if `x` is a list).
#' @param kaiser_line Logical. If `TRUE`plots a line at y = 1.
#' @param ylab Character. Y-axis label.
#' @param addtitle Logical. If `TRUE`, adds a plot title.
#' @param prn Logical. If `TRUE`, prints the plot.
#' @param rtn Logical. If `TRUE`, returns the ggplot object.
#'
#' @return If `rtn = TRUE`, a `ggplot` object; otherwise `NULL` (invisibly).
#'
#' @export
 screeplot = function(x, nplot, kaiser = F, ylab = "eigenvalues", addtitle = T, prn = TRUE, rtn = FALSE){
  
  if(missing(nplot)) nplot = length(x)
  if(is.list(x)){
    if(any(names(x) == "values"))
      x = x$values
    else
      stop("x must be a list with elment `$values', or a vector of eigenvalues")
  }
  df = data.frame(order = 1:nplot,
                  eigenvalue = x[1:nplot])
  scree_pl = ggplot(df, aes(x = order, y = eigenvalue)) + geom_point(size = 2) + geom_line() + labs(y = ylab) + theme_tuto()
  if(kaiser_line) scree_pl + ggplot2::abline(h = 1)
  
  if(addtitle) scree_pl = scree_pl + labs(title = "screeplot") +
    theme(plot.title = element_text(hjust = 0.5))
  
  if(prn)
    print(scree_pl)
  if(rtn)
    return(scree_pl)
}

#not distributed coz needs cowplot
#both scree and wach==========
#scree_watc_plot = function(eigvals, nplot, p, n, gamma, wa_cor = T, wa_kai = F, screeperc = TRUE, scree_ylab = "variance explained", wa_nfit_line = NULL){
#   if(missing(p)) p = length(eigvals)
#   if(missing(n))
#     stop("must pass n in scree_watc_plot")
#   
#   if(sum(eigvals) == p) wa_cor = T
#   sc = screeplot(eigvals, nplot = nplot, perc = screeperc, ylab = scree_ylab, rtn = T, prn = F)
#   wc = wachterqq(eigvals = eigvals, p = p, n = n, cor = wa_cor, nplot = nplot, nfit_line = wa_nfit_line, addtitle = TRUE, prn = FALSE, rtn = TRUE)
#   cowplot::plot_grid(sc, wc)
# }

# Rounds  a list
#
# prints or return a rounded list
#
# @param li a list of numerical objects
# @param d nuber of digits
# @return rounded list
roundl = function(li, d = 2) lapply(li, round, digits = d)


# Change the sign of selected sparse components
#
# Multiplies by \eqn{-1} the loadings, contributions, and (when present) scores and
# component-correlation entries for the components listed in `index_to_change`.
# This is useful because (sparse) principal components are defined up to sign.
#
# @param spca_obj An object of class `spca`.
# @param index_to_change Integer vector of component indices whose sign should be flipped.
#
# @return The modified `spca_obj`, with the selected components sign-flipped.
#
# @noRd
# @export
changeSign_loads_spca = function(spca_obj, index_to_change){
  
  n = length(index_to_change)
  for (i in 1:n){
    spca_obj$loadings[, index_to_change[i]] = - spca_obj$loadings[, index_to_change[i]]
    spca_obj$contributions[, index_to_change[i]] = - spca_obj$contributions[, index_to_change[i]]
    
    if(!is.null(spca_obj$loadingslist)){
      spca_obj$loadingslist[[index_to_change[i]]] =
        -spca_obj$loadingslist[[index_to_change[i]]]
    }
    if (!is.null(spca_obj$scores))
      spca_obj$scores[, index_to_change[i]] = - spca_obj$scores[, index_to_change[i]]
    if (!is.null(spca_obj$corComp)){
      spca_obj$corComp[index_to_change[i], ] = -spca_obj$corComp[index_to_change[i], ]
      spca_obj$corComp[, index_to_change[i]] = -spca_obj$corComp[, index_to_change[i]]
    }
  }
  return(spca_obj)
}



# makes VIF for an spca object
#
# experimental use with care, uses only vifSC because vifXC  slow
#
# @param spca_obj surprisingly, an spca object
# @param M data matrix or variance/correlation matrix (same result)
# @param intercept if include intercept, use TRUE only for data matrix not centered.
#
# @return a list of vectors with vif values for each component
make_vif.spca = function(spca_obj, M, intercept = FALSE){
  stopifnot(any(class(spca_obj) == "spca"))
  if(is.null(colnames(M))){
    namx = paste("Var", 1:ncol(M))
  }
  else{
    namx = colnames(M)
  }
  viff = as.list(1:spca_obj$ncomps)
  if(!isSymmetric(M))
    M = cor(M)
  for (j in 1:spca_obj$ncomps){
    if(spca_obj$cardinality[j] > 1)
      viff[[j]] = vifSC(M[spca_obj$ind[[j]], spca_obj$ind[[j]]], 1:spca_obj$cardinality[j])
  }
  for (i in 1:spca_obj$ncomps){
    names(viff[[i]]) = namx[spca_obj$ind[[i]]]
    ##   viff[[i]] = sort(viff[[i]])
  }
  return(viff)
}

makevec = function(vec, n){
  le = length(vec)
  if(le == 0)
    stop("need to pass a vector of lenght > 0")
  if(le < n){
    m = n - le
    v = c(vec, rep(vec[le], m))
  }
  if (le == n)
    v = vec
  return(v)
}

anyna = function(x) any(is.na(x))