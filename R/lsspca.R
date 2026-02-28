

# tutorial 2025
## requires C++ func for matrix comp =================P

#lsspca=====================
#' \strong{Computes LS SPCA using different selection methods}
#' 
#' For each component, the variables are selected so as to explain
#' a percentage \emph{alpha} of the vexp by the corresponding principal component.
#' @param X dataframe or matrix with data. 
#' @param alpha = 0.95. Real in [0,1] (optional if ncbyvexp < 1). Minimum \code{R^2} of the regression of the residual PC for variable selection. 
#' @param ncomps = 0 Number of components to compute (if passed 0, reset to the number of variables)
#' @param ncbyvexp = 1 (optional if ncomps < 1) real in [0, 1] compute components until this fraction of total vexp explained. whichever satisfied first ncomps.
#' @param method a character vector, one of "p", "c", "u", specifying how the lsspca loadings are computed for each component. The last value will be applied to all remaining components.
#' @param  varselection one of "stepwise", "backward" or "forward". Default "stepwise". Initial letter enough.
#' @param maxcard integer (optional) a vector or one integer. Missing values filled with last value.
#' @param force_in NULL or list or vector of indices must be in component. 
#' @param force_out NULL or list or vector of indices cannot be in component.  
#' @param scalex = FALSE if TRUE the variables are scaled to unit length. Variables are automatically scaled to zero mean (with warning) if they are not.
#' @param mkvif = FALSE, if true computes the variance inflation for each variable selected  in every component May slowdown.
#' @details
#' \code{alpha} controls the $R^2$. The function may be slow because it uses the package leaps for variable selection. This is to keep results the same as in the paper.  
#' 
#' @return an \code{\link{spca_object}}.
#' @export
lsspca =
  function (X, alpha = 0.95, ncomps = 0, ncbyvexp = 1, method = "p", varselection = c("stepwise", "backward", "forward"), maxcard = 0, force_in = NULL, force_out = NULL, scalex = FALSE, mkvif = FALSE) 
  {
  
  if(is.data.frame(X)) X = as.matrix(X)
  p = ncol(X)
  n = nrow(X)

  
## validation ============================
  
    
    if (ncomps == 0) {
      if ((ncbyvexp == 1)) 
        (stop("need give either ncomps > 0 or ncbyvexp , 1"))
      else ncomps = p
    }
    else {
      ncbyvexp = 1
    }
  
# method ------------------
  
  for(j in 1:length(method)){
  ssearch = substr(method[j], 1, 1)
  method[j] = switch(ssearch,
                        p = "p",
                        c = "c",
                        u = "u",
                        stop("method must be one of p, c, or u")
  )
  }
  # fills method up to p
  if (length(method) < p)
    method = makevec(method, p)
  
    if (is.vector(maxcard)) {
      maxcard[maxcard == 0] = p
      if (length(maxcard) < ncomps) 
        lm = length(maxcard)
        maxcard = c(maxcard, rep(p, ncomps - lm))
      }
    else {
      stop("must pass a vector or an integer as maxcard")
    }
    ssearch = substr(varselection[1], 1, 1)
    varselection = switch(ssearch,
                        b = "backward",
                        f = "forward",
                        s = "seqrep",
                        stop("varselection must be one of stepwise, backward or forward")
    )


# converts force in/out to lists =========================
# at the end check content of of each to avoid weird   
    
    if(!is.list(force_out))
      if(is.vector(force_out))
        force_out = list(force_out)
    
    a = vector("list", ncomps); 
    for(i in 1:p) a[i] = force_out[i]#unlist(force_out[i]) 
    force_out = a
    
    if(!is.list(force_in))
      if(is.vector(force_in))
        force_in = list(force_in)
    
    a = vector("list", ncomps); 
    for(i in 1:p) a[i] = force_in[i]#unlist(force_in[i]) 
    force_in = a
    
    chk = sapply(1:ncomps, function(i, x, y) length(intersect(x[[i]], y[[i]])) > 0, x = force_in, y = force_out)
    if(any(chk)){
      stop(paste("cannot force in and out same variable. check elements", which(chk)))
      }

    if (is.null(colnames(X))) 
      namx = paste("V", 1:p)
    else namx = colnames(X)
    if (scalex == TRUE) 
      X = scaleC(X)
    else if (any(abs(colMeans(X)) > 10^-6)) {
      message("Centering the columns of X to zero mean")
      X = scaleC(X, scale = FALSE)
    }

#  working objects ===========    
    card = rep(0, ncomps)
    ind = vector("list", ncomps)
    vexp = rep(0, ncomps)
    cvexp = rep(0, ncomps)
    A = matrix(0, p, ifelse(ncomps == 0, p, ncomps))
    colnames(A) = paste0("PC", 1:ncomps)
    contributions = A
    loadlist = vector("list", ncomps)
    scores = matrix(0, n, ncomps)
    colnames(scores) = paste0("PC", 1:ncomps)
    
    if (any(method == "u"))
      R = matrix(0, p, ncomps)
    
    K = X
    nc = 0
    j = 1
    stopComp = FALSE
    
# looping ===================   

  ## till r2_min or ncombyvexp satisfied    
    while (stopComp == FALSE) {
#     print(paste("start comp  ", j))
      
      D = ataC(K)
      r_ee = EigenC(D)
      pc = avC(K, r_ee$vec[, 1, drop = T])
    #  pc = K %*% r_ee$vec[, 1]
      if (j == 1) {
        Cmat = D
        vexpPC = r_ee$val
        vtot = sum(r_ee$val)
  }

# Xd with variable selected ==================    
#   reduce matrix and indices if force_out[[j]] not null 
  if (length(force_out[[j]]) > 0) {
    indall = setdiff(seq_len(p), force_out[[j]])
    Xd = X[, indall, drop = FALSE]
    } 
  else {
    indall = seq_len(p)
    Xd = X
  }
  
  # Map force_in (given in original X indices) to Xd indices
      force_in_j = force_in[[j]]
      if (length(force_in_j) > 0) {
        force_in_Xd = match(force_in_j, indall)
        if (anyNA(force_in_Xd)) {
          stop(sprintf("force_in[[%d]] contains indices removed by force_out[[%d]]", j, j))
        }
      } 
      else {
        force_in_Xd = NULL
      }    
      
      ##Reg search=================  
      # leaps search
      ssr = leaps::regsubsets(x = Xd, y = pc, method = varselection[1], nbest = 1, force.in = force_in_Xd, nvmax = maxcard[j], intercept = FALSE
      )# THIS is ok coz pc from K
      
      aa = summary(ssr)
      mrsq = any(aa$rsq >= alpha)  
      if(mrsq){
        if(!is.null(ssr$force_in)) { #leaps moves forced in variables at the beginning
          w  = aa$which[indmodel, ]
          
          # Check: in regsubsets' own naming, force.in must be TRUE
          stopifnot(all(w[ssr$force.in]))  # ssr$force.in is stored in the same ordering    
          # Convert selected columns back to original X indices via names
          sel_Xd   = match(names(w)[w], colnames(Xd))     # positions in Xd
          indmodel    = indall[sel_Xd]    
          
        }
        else{
          if (method[j] == "u") 
            indmodel = which((aa$rsq >= alpha) & ((1:length(aa$rsq)) >= j))[1]
          else
            indmodel = which(aa$rsq > alpha)[1]
          
          ind[[j]] = indall[aa$which[indmodel, ]]
          card[j] = length(ind[[j]])
        }
      }
        else{
          warning(paste("Component", j, "could not reach ", alpha, "VEXP\n", "including all variables"))
          ind[[j]] = indall
        }# end leaps
# Compute coefficients ====================  
      
      Xd = X[, ind[[j]]]
        uspcafail = FALSE
#Only need compute loadings if more than one variable selected        
        if (length(ind[[j]]) > 1) {
          Sd = ataC(Xd)
          if (method[j] == "u") {
            if (j == 1) {
              M = aatC(crossprod(Xd, X))
              ga = GenEigenC(M, Sd)
              a = ga$vec[, 1]
            }#end i = 1
            else {#j >1
              # if singular skips and go to cspcacspca         
              if(!rank_sym(Sd)){
                warning(paste("pspca: Comp", j, "Cannot compute uncorrelated component, switch to uncorelated"))
                uspcafail = TRUE
              }#end chck rank
              
              else{
                xxd = atbC(X, Xd)
                H = atbC(A[, 1:(j - 1), drop = FALSE], xxd)
                Sm = solveC(Sd)
# this was for W = abtc(E, H)    E = H %*% Sm
                W = abtC(abC(H, Sm), H)
                # if singular break and go to cspca         
                if(!rank_sym(W)){
                  warning(paste("pspca: Comp", j, "Cannot compute uncorrelated component, switch to uncorelated"))
                  uspcafail = TRUE
                  #break
                }#end chck rank
                else{
                  G = diag(card[j]) - atbC(H, solveC(W)) %*% H %*% Sm
                  M = ataC(abtC(xxd, G))
                  ga = GenEigenC(M, Sd)
                  a = ga$vec[, 1]
                }#end else chck prod mat singular
              }#end else chck Sd singular
            }#end "u"
            }#end method[j] = u
          if((method[j] != "u") | uspcafail) {#this don't need
# CSPCA            
            if ((method[j] == "c") | uspcafail) {
              M = ataC(atbC(K, Xd))
              ga = GenEigenC(M, Sd)
              a = ga$vec[, 1]
            }#end method[j] = c
# PSPCA already done in search, actually
          if ((method[j] == "p")) {
              alm = lm(pc ~ Xd - 1)
              a = alm$coefficients
            }#end method[j] p
          }
        }#end compute loadings 

## compute stats and stores step j results=========          
      if(length(a) > 1){   
       a = a/sqrt(sum(a^2))
          if (all(a <= 0)) 
            a = -a
          scores[, j] = Xd %*% a
          bb = stats::cor(pc, scores[, j])
          if (bb < 0) {
            a = -a
            scores[, j] = -scores[, j]
          }
        }#end length(a)> 1
        else {
          a = 1
          scores[, j] = X[, ind[[j]]]
          bb = stats::cor(pc, scores[, j])
          if (bb < 0) {
            a = -a
            scores[, j] = -scores[, j]
          }
        }#end else length(a) >1 
      names(a) = namx[ind[[j]]]
      A[ind[[j]], j] = a
      contributions[ind[[j]], j] =  a/sum(abs(a))
      loadlist[[j]] = a
 ## Deflation and stop =====================================      
      tmp = deflXCforR(c(scores[, j]), K)
      vexp[j] = tmp$vexp
      cvexp[j] = ifelse(j == 1, vexp[j], sum(vexp[1:j]))
      K = tmp$K
      nc = nc + 1
      #print(paste("done comp", j))
      if ((cvexp[j] >= ncbyvexp * vtot) | (j == ncomps) | (j == 
          p)) {
        stopComp = TRUE
        }
      else {
        j = j + 1
      }
    }#end while
    
    # if(nc != ncomps){
    #   print(nc)}
#    browser()
# output ==================   
    ncomps = nc
    A = A[, 1:ncomps]
    aa = makeVexpC(A, Cmat) 
    vexp = aa$vexp
    cvexp = aa$cvexp
    rownames(contributions) = namx
    rownames(A) = namx
    ind = lapply(ind, function(x, na) {names(x) = na[x]; x}, na = namx) 
    out = list(loadings = A, 
               contributions = contributions[, 1:ncomps], 
               vexp = vexp[1:ncomps]/vtot, 
               vexpPC = vexpPC[1:ncomps]/vtot, 
               cvexp = cvexp[1:ncomps]/vtot, 
               rcvexp = cvexp[1:ncomps]/cumsum(vexpPC[1:ncomps]),
               rpcvexp =  vexp[1:ncomps]/vexpPC[1:ncomps],
               scores = scores[, 1:ncomps], 
               ncomps = ncomps, 
               ind = ind[1:ncomps], 
               cardinality = card[1:ncomps], 
               loadlist = loadlist[1:ncomps], 
               method = method
               )
    if (ncomps > 1) {
      if (any(method[j] != "u")) {
        out$corComp = MkCorCompMat(A[, 1:ncomps], stats::cor(X), d = ncomps)
      }
    }
    out$Call = match.call()
    class(out) = c("spca", "list")
    if(mkvif){
      out$vif = make_vif_R(out, X, prn = FALSE)
        for (i in 1:ncomps) {
        names(out$vif[[i]]) = namx[ind[[i]]]
        out$vif[[i]] = sort(out$vif[[i]])
      }
    }
    return(out)
  }

rank_sym = function(A, tol = NULL){
  ev = eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if (is.null(tol)) tol = max(dim(A)) * .Machine$double.eps * max(abs(ev))
  all(abs(ev) > tol)
}