

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

#x = hsp_scale_list
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

#rm(aggregate_by_scale)

# aggregate_by_scale = function(x, ind, only.nonzero = T, addScaleNames = T) {
#   if(is.vector(x))
#     out = tapply(x, ind, sum)
#   else{
#     out = apply(x, 2, function(y, ii) tapply(y, ii, sum), ii = ind)
#   if(only.nonzero)
#     out = drop.levels(out[rowSums(abs(out)) > 0, ])
#   }
#   return(out)
# }

##scale_y_continuous(labels = scales::percent)

make_colours = function(n, pal = c("cbb", "ggplot")){
  if (grepl("^c", pal[1])){
    cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    if (n <= 8) 
      cols = cbb[1:n]
    else{
      warnings("cannot use bcc pallette with more than 8 colours. Switch to ggplot")
      pal = "ggplot"
    }
  }
  if(grepl("^p", pal[1])){
    if (n < 10){
      pf.pal <- rev(c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
      cols = pf.pal[1:n]
    }
    else{
      warnings("cannot use pf pallette with more than 8 colours. Switch to ggplot")
      pal = "ggplot"
    }
  }
  if (grepl("^g", pal[1])){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols <- gg_color_hue(n)
  }
  return(cols)
}
# make_colours(7, "p")
# make_colours(7, "cbb")
# make_colours(7, "gg")

# pl <- pl + ggplot2::scale_fill_manual(values = make_colours(n))

mkVexpTuto25 = function(A, S, totv, cvex = T){
  
  p = ncol(A)
  M = S %*% A
  cvexp = rep(0, p)
  cvexp[1] = sum(M[, 1]^2)/(crossprod(A[, 1], M[, 1]))
  if(p > 1){
    for(i in 2:p){
      cvexp[i] = sum(diag(M[, 1:i] %*% solve(crossprod(A[, 1:i], M[, 1:i])) %*%  t(M[, 1:i])))
    }
  }
  cvexp = cvexp/totv
  vexp = c(cvexp[1], diff(cvexp))
  if(cvex)
    return(list(vexp = vexp, cvexp = cvexp))
  else
    return(vexp = vexp)
}
#wachter=============
#scale_y_continuous(labels = scales::percent)

wachterqq = function(eigvals, p, n, gamma, cor = T, nplot, nfit_line = NULL, addtitle = TRUE, prn = TRUE, rtn = FALSE){
  
  if(missing(gamma)) gamma = n/p
  
  if(missing(nplot)) nplot = length(eigvals)
  
  probs <- ((p - (1:p) + 1)- 0.5)/p
  mp_quantiles <- RMTstat::qmp(p = probs, svr = gamma)
  
  if(cor) mp_quantiles <- p * mp_quantiles/sum(mp_quantiles)
  
  df = data.frame(expected = mp_quantiles[1:nplot], observed = eigvals[1:nplot])
  pl = ggplot(df, aes(x = expected, y = observed)) + geom_point(size = 2)
  
  if ((is.numeric(nfit_line)) && (nfit_line > 0))
    pl = pl + geom_smooth(data = df[(nplot - nfit_line):nplot, ], se = F, method = "lm")
  if(addtitle)
    pl = pl + labs(title = "wachter qq-plot")+
    theme(plot.title = element_text(hjust = 0.5))
  
  if(prn)
    print(pl)
  if(rtn)
    return(pl)
}

#myscreeplot =====================
myscreeplot = function(x, nplot, perc = TRUE, ylab = "variance explained", addtitle = T, prn = TRUE, rtn = FALSE){
  
  if(missing(nplot)) nplot = length(x)
  if(is.list(x)){
    if(any(names(x) == "values"))
      x = x$values
    else
      stop("x must be a list with elment `$values', or a vector of eigenvalues")
  }
  if(perc) x = x/sum(x)
  df = data.frame(order = 1:nplot,
                  vexp = x[1:nplot])
  scree_pl = ggplot(df, aes(x = order, y = vexp)) + geom_point(size = 2) + geom_line() + labs(y = ylab)
  if(addtitle) scree_pl = scree_pl + labs(title = "screeplot") +
    theme(plot.title = element_text(hjust = 0.5))
  
  if (perc)
    scree_pl = scree_pl + scale_y_continuous(labels = scales::percent)
  
  if(prn)
    print(scree_pl)
  if(rtn)
    return(scree_pl)
}
#both scree and wach==========
scree_watc_plot = function(eigvals, nplot, p, n, gamma, wa_cor = T, wa_kai = F, screeperc = TRUE, scree_ylab = "variance explained", wa_nfit_line = NULL){
  if(missing(p)) p = length(eigvals)
  if(missing(n))
    stop("must pass n in scree_watc_plot")
  
  if(sum(eigvals) == p) wa_cor = T
  sc = myscreeplot(eigvals, nplot = nplot, perc = screeperc, ylab = scree_ylab, rtn = T, prn = F)
  wc = wachterqq(eigvals = eigvals, p = p, n = n, cor = wa_cor, nplot = nplot, nfit_line = wa_nfit_line, addtitle = TRUE, prn = FALSE, rtn = TRUE)
  cowplot::plot_grid(sc, wc)
}

IMPROVE DOCUMENTATION

#' @param smpc either spca object or matrix or data.frame of contributions
#' @param vargroups list or vector or factor of grouping indices
#' @param ncomp how many components to consider, default all 

#make_contByscale==========================
make_contByscale = function(smpc, vargroups, ncomp, only.nonzero = T){
  if(is.list(vargroups)) 
    vargroups = unlist(vargroups)
  if(is.character(vargroups))
    vargroups = forcats::as_factor(vargroups)
  if(any(class(smpc) == "spca")) 
    cont = smpc$contributions
  else
    if (is.matrix(smpc) || is.data.frame(smpc))
      cont = as.matrix(smpc)
    else
      stop("smpc must be either spca object or matrix or dataframe" )
    
    if(missing(ncomp))    ncomp = ncol(cont)
    
    ContByscale = sapply(1:ncomp, function(i, l, s) tapply((l[, i]), s, sum), l = cont, s = vargroups)
    
    colnames(ContByscale) = paste0("Comp", 1:ncomp)
    
    if(only.nonzero)
      ContByscale = ContByscale[rowSums(ContByscale) != 0, ]
    return(ContByscale)
}



#' Rounds  a list
#'
#' prints or return a rounded list
#'
#' @param li a list of numerical objects
#' @param d nuber of digits
#' @return rounded list
#' @export
roundl = function(li, d = 2) lapply(li, round, digits = d)



#' changes the signs to loadings in spca object
#'
#' laodings, contributions, correlation between components and scores are changed
#'
#' @param spca_obj an object of class spca
#' @param index_to_change vector of indices of loaidngs to change sign
#' @return  the spca object with changed signs
#' @export
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



#' makes VIF for an spca object
#'
#' experimental use with care, uses only vifSC because vifXC  slow
#'
#' @param spca_obj surprisingly, an spca object
#' @param M data matrix or variance/correlation matrix (same result)
#' @param intercept if include intercept, use TRUE only for data matrix not centered.
#'
#' @return a list of vectors with vif values for each component
#' @export
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
      viff[[j]] = spcaTutoPack::vifSC(M[spca_obj$ind[[j]], spca_obj$ind[[j]]], 1:spca_obj$cardinality[j])
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