## ok from devel
makevexp =  function (A,D) {
{ ## should add a flag for unc 
  ## new improved using ind
  ## revised to allow one component only
  ## computes vexp for uncorrelated components
}
if (is.matrix(A)){
  A =  t(t(A)/ sqrt(colSums(A^2)))
  nd = ncol(A)
  vt = sum(diag(D))
  vexp = rep(0, nd)
  for (i in 1:nd) {
    ind = which(abs(A[,i])> 0.001)
    if (length(ind) == 1){
      B = D[,ind] 
      b = D[ind,ind]
      vexp[i] = sum(B^2)/(vt * b)
    }
    else{
      B = D[,ind] %*%  A[ind, i]
      b = drop(crossprod(B[ind,], A[ind, i]))
      vexp[i] = crossprod(B, B)/(vt * b)
    }
    D = D - tcrossprod(B ,B)/b
  }
}
else{
  ind = which(abs(A) > 0.001)
  A = A[ind]/sqrt(sum(A[ind]^2))
  vt = sum(diag(D))
  B = D[,ind] %*% as.matrix(A)
  b = drop(crossprod(B[ind,], as.matrix(A)))
  vexp = crossprod(B, B)/(vt * b)
}
return(vexp)
}

eig = function(D, nd = 4, prn = FALSE){
  ## does eigendecomposition returning also vexp and cvexp
  ## returns object spca  
  ee = eigen(D, symmetric = TRUE)
  ee$vexp = ee$val/sum(ee$val)
  ee$cvexp = cumsum(ee$val)/sum(ee$val)
  if (prn == TRUE){
    out = list(A = ee$vec[,1:nd], vexpv = ee$vexp[1:nd])#, vexp = ee$vexp[1:nd])
    class(out) = "spca"
    print(out)
  }
  return(ee)  
}

make.uncLoad = function(A, S){
  ## ortogonalises loadings
  p = ncol(S)
  if (is.vector(A))
    d = 1
  else
    d = ncol(A)
  if (d > 1){
    B = A   
    Z = diag(1, p)
    for (i in 2:d){
      Z = makez(B[,i-1], S, Z)
      B[,i] = Z %*% A[,i]
    }  
    B = t(t(B)/ sqrt(colSums(B^2)))
    return(B)
  }
  else
    return(A)
}

make.corx = function(S, A){
  ## computes the cor between the components and each x-variable  
  M = S %*% A
  top = t(t(M)/diag(sqrt(crossprod(A,M))))
  top = top/sqrt(diag(S))
  return(top)  
} 

make.cor = function(D, A, dgt = 4){
  ## computes correlation among components
  if (is.vector(A))
    out = (NULL)
  else 
    if (ncol(A) == 1)
      out = (NULL)
  else 
    if (is.matrix(A) & ncol(A) > 1 ){  
      out = crossprod(A, D) %*% A
      o = sqrt(diag(out))^-1
      out = round(diag(o) %*% out %*% diag(o), dgt)
    }
  return(out)
}

makez = function(a, S, Z){
  ### updates the Z matrix with next loadings vector
  ## only for unc components
  p = ncol(S)
  if (!is.vector(a))
    if(ncol(a)> 1)
      stop("you must pass a vector or column matrix to makez, consider using makezM")
  a = as.matrix(a)
  if (missing(Z))
    diag(p) - (a %*% crossprod(a,S))/(drop( crossprod(a,S) %*% a)) 
  else  
    Z - (a %*% crossprod(a,S))/(drop( crossprod(a,S) %*% a))  
}

makezM = function(A,D){
  ## makes Z matrix from a full set of loadings, 
  ## must use this with correlated components
  diag(ncol(D)) - A %*% solve(crossprod(A,D) %*% A) %*% crossprod(A,D) 
}

make.vexpeig = function(ee, ndim, cum = TRUE, dgt = 3){
  ## computes the proportional variance explained by eigenvectors
  if(missing(ndim))
    ndim = length(ee$val)
  if (cum == TRUE)
    out = round(cumsum(ee$val[1:ndim])/sum(ee$val), dgt)
  else
    out = round((ee$val[1:ndim])/sum(ee$val), dgt)
  return(out)  
}

make.cont = function(smpc){
  ## standardise a matrix of loadings to unity l1 norm
  ## scales loadings to unit L1 contributions
  ## input either a spca object or matrix of loadings
  ## v2 modified, doesnt return vexp anymore
  if (any(class(smpc) == "spca")){
    if (any(names(smpc) == "contributions"))
      Ac = smpc$contributions
    else{
      Ac = t(t(smpc$loadings)/colSums(abs(smpc$loadings)))
    }
  }  
  else{ 
    if (is.list(smpc) & length(smpc) == 1)
      smpc = smpc[[1]]
    if ( is.matrix(smpc))
      Ac = t(t(smpc)/colSums(abs(smpc)))

    else
      if (is.vector(smpc))
        Ac = smpc/sum(abs(smpc))
    else
      stop("a matrix of loadings or an spca object is needed")      
  }
  return(Ac)
}

get.minload = function(smpc, perc = FALSE, eps = 0.001){
  ## returns the non-zero loading with the smallest absolute value for each column
  ## input a matrix of loadings or an spca object
  ## if perc == TRUE will compute on the percent contribution  
  if (any(class(smpc) == "spca"))
    smpc = smpc$loadings
  else 
    if (! is.matrix(smpc))
      stop("get.minload: a matrix of loadings or an spca object is needed")  
  if (perc == TRUE)
    smpc = make.cont(smpc)
  d = ncol(smpc)
  gl = function(x)
    min(abs(x[abs(x)> eps]))
  apply(smpc, 2, gl)
}

get.card = function(A, eps = 0.01){
  ## returns the cardinality of the columns of a matrix of loadings  
  if (any(class(A) == "spca"))
    A = A$loadings
  if(is.null(dim(A)))
    sum(abs(A) > eps)
  else
    colSums(abs(A)> eps)
}

mult.eigen = function(A, b, ind, power = 1L){
{#multiplies A %*% diag(b) %*% t(A)
  # used for fitted variance matrix
}
if (missing(ind))
  ind = 1:ncol(A)
if (length(b[ind]) != ncol(A[,ind]))
  stop("matrix A and vector b must be of compatible dimension")
if (power == 1)
  t(A[,ind] %*% (t(A[,ind]) * b[ind]))
if (power != 1)
  if (power == 0.5 & any(b < 0))
    stop("mult.eigen: cannot take sqrt of negative b")
t(A[,ind] %*% (t(A[,ind]) * (b[ind])^power))
}

myprintspca = function(smpc, cols, digits = 3, rows, noprint = 1E-03, 
                       rtn = FALSE, perc = FALSE, namescomp = NULL){
  ### generic print for spca objects 
  ## used by spca.comp
  ## v2.1 fixed problem with only 1 component
  if (any(class(smpc) == "spca")){  
    if (perc == TRUE){
      if(any(names(smpc) == "contributions"))
        A = smpc$contributions
      else
        A = make.cont(smpc$loadings)
    }
    else{
      if (is.null(smpc$As))
        A = smpc$loadings
      else 
        A = smpc$As
    }
  }
  else{#
    if (is.matrix(smpc) | is.vector(smpc)){
      A = make.cont(smpc)
    }
    else
      stop(paste("The argument must be either a smpc object or a matrix, not a", class(smpc)))
  }
  
  if (is.vector(A))
    A = as.matrix(A)
  ## end if
  if (missing(cols))
    cols = 1:ncol(A)
  else
    if (length(cols) == 1L)
      cols = 1:cols
  
  if (missing(rows))
    rows = 1:nrow(A)
  else
    if (length(rows) == 1L)
      rows = 1:rows
  
  A = as.matrix(A[rows,cols])
  ## assigns names to laodings
  if (!is.null(namescomp) & length(namescomp) == ncol(A)){
    colnames(A) = namescomp
  }
  else{
    if (!is.null(namescomp) & length(namescomp) != ncol(A))
      warning("the length of namescomp is incorrect, automatic names assigned")
    colnames(A) = paste("Comp",1: ncol(A), sep = "")
  }
  
  
  
  if (perc == TRUE)
    fx <- format(round(A*100, max(digits-2,0)),drop0trailing = TRUE, justify = "centre")
  else
    fx <- format(round(A, digits),drop0trailing = TRUE, justify = "centre")
  names(fx) <- NULL
  nc <- nchar(fx[1L], type = "c")
  fx[abs(A) < noprint] <- paste(rep(" ", nc), collapse = "")
  #  ind = (abs(A)> noprint & abs(A) < noprint)
  #  fx[ind] = "--"
  fx = format(fx, justify = "right" )
  if (any(class(smpc) == "spca")){
    ## chv
    vexp = smpc$vexp[cols]
    doo = rep("-----", ifelse(is.null(ncol(fx)), 1, ncol(fx)))
    fx = rbind(fx, doo, round(100*vexp,1))
    rownames(fx)[nrow(fx)-1] = ""    
    rownames(fx)[nrow(fx)] = "PCVE"
  }  
  if (perc == TRUE)
    print("Percentage Contributions")
  else
    print("Loadings")
  if (ncol(A) == 1L){
    #    fx = t(fx)
    print(t(fx), quote = FALSE)#, ...)
    
    
    #     rownames(fx)[1] = "Loadings"
  }
  else
    print(fx, quote = FALSE)#, ...)
  
  if(rtn == TRUE){   
    return(fx)    
  }  
  else 
    invisible()
}


## computes variance inflation factors for all components
# \strong{Computes variance inflaction factors for variables in spca components.}
# 
# result provided by SPCA
# 
# @param spca_obj An spca object, typically from SPCA.
# @param M the data or covariance/correlation matrix
# @param intercept Logical: should the intercept be included in the model?
# @param pseudo, should compute vifs using generalised inverse?
# Only when M is data matrix.
# @param prn Logical. Should the function print the results?
# @param digits. Number of decimal figures to print and return. See details.
# @details A wrapper for vifSC in C++. May run into problems if the variables
# are perfectly multicollinear, in that case use \code{pseudo = TRUE}.
# Vifs values are named with \code{colnames(M)}.
# If \code{digits} is 0, loadings are returned unrounded and printed (if enabled)
# to two digits. If \code{digits} > 0, loadings are returned rounded to that value.
# @return A list of ordered vifs for each component.
# # was method for class. no nee

make_vif_R = function(spca_obj, M, intercept = FALSE, pseudo = FALSE, 
                         prn = TRUE, digits = 0){
  stopifnot(any(class(spca_obj) == "spca"))
  if(is.null(colnames(M))){
    namx = paste("Var", 1:ncol(M))
  }
  else{
    namx = colnames(M)
  }
  if(is.null(spca_obj$ncomps)){
    spca_obj$ncomps = min(length(spca_obj$vexp), ncol(spca_obj$loadings))
    warning(paste("ncomps set to ", spca_obj$ncomps))
  }
  if(is.null(spca_obj$cardinality)){
    spca_obj$cardinality = colSums(abs(spca_obj$loadings) > 0.01)
    warning(paste("cardinality set to  ", spca_obj$cardinality))
  }
  if(is.null(spca_obj$ind)){
    spca_obj$ind = lapply(1:ncomps, function(i, A, eps = 0.01) which(abs(A[, i]) > 0.01), 
                          A = spca_obj$loadings)
  }
  viff = as.list(rep(NA, spca_obj$ncomps))
  if(!isSymmetric(M))
    M = cor(M)
  for (j in 1:spca_obj$ncomps){
    if(spca_obj$cardinality[j] > 1)
      if(pseudo == FALSE)
        viff[[j]] = vifSC(M[spca_obj$ind[[j]], spca_obj$ind[[j]]], 1:spca_obj$cardinality[j])
    else
      viff[[j]] = vifSPseudoC(M[spca_obj$ind[[j]], spca_obj$ind[[j]]], 1:spca_obj$cardinality[j])
  }
  for (i in 1:spca_obj$ncomps){
    names(viff[[i]]) = namx[spca_obj$ind[[i]]]
    viff[[i]] = sort(viff[[i]])
    if (digits > 0)
      viff[[i]] = round(viff[[i]], digits)
  }
  if(prn == TRUE){
    for (i in 1:length(viff)){
      cat(paste("Comp", i, "\n"))
      print(round(viff[[i]], max(2, digits)))
      cat("\n")
    }
  }
  return(viff)
}
