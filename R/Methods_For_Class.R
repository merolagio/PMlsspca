
#is.pca==============
#' \strong{Verifies if an object is of class spca}
#' 
#' @param x Any object suspected of being of class spca.
#' @details The function carries out checks for the minimal components of an
#' \emph{spca} object. It checks the presence and mode of
#'  a matrix of loadings named \emph{loadings} 
#'  a vector of variance explained by the components named \emph{vexp} 
#'  a vector of variance explained by the PCss named \emph{vexpPC}
#' an integer declaring the number of components named \emph{ncomps}.
#' If any of these is missing the object is deemed as not \emph{spca} and a
#' warning is issued, even if the object's class is \emph{spca}.
#' @return Logical: TRUE if object is of class spca and contains the minimal configuration, 
#' FALSE otherwise. 
#' @export 
is.spca <- function(x){
  isclass = inherits(x, "spca")
  isloadings = (!is.null(x$loadings))
  if (isloadings)
    isloadingsmat = is.matrix(x$loadings)
  
  isvexp = (!is.null(x$vexp))
  if (isvexp)
    isvexpvec = is.vector(x$vexp)
  
  isvexpPC = (!is.null(x$vexpPC))
  if (isvexpPC)
    isvexpPCvec = is.vector(x$vexpPC)
  
  isncomps = (!is.null(x$ncomps))
  if (isncomps)
    isncompsvec = is.vector(x$ncomps)
  
  isit = TRUE
  out = list()
  if (isclass){
    out = c(out, "the object is of class spca")
  }
  else{ 
    out = c(out, "the object is not of class spca")
    isit = FALSE
  }
  if (!(isloadings && isloadingsmat) ){
    out = c(out, "it must contain a matrix of loadings named loadings")
    isit = FALSE
  }
  
  if (!(isvexp && isvexpvec) ){
    out = c(out, "it must contain a vector of variances explained by the components named vexp")
    isit = FALSE  
  }
  
  if (!(isvexpPC && isvexpPCvec) ){
    out = c(out, "it must contain a vector of variances explained by the PCs named vexpPC")
    isit = FALSE  
  }
  
  if (!(isncomps && isncompsvec) ){
    out = c(out, "it must contain an integer declaring the number of components named ncomps")
    isit = FALSE  
  }
  if(!isit){
    for (o in out)
      warning(o)
  }
  return(isit)
}
# print.spca =======================
#' \strong{Prints the sparse loadings from an spca object}
#' 
#' Prints sparse loadings omitting the zero ones and giving the cumulative
#' variance explained.
#' 
#' 
#' @param spca_obj An spca object.
#' @param cols A vector indicating which components should be printed. Default
#' all. If an iteger is passed, it is set to 1:cols.
#' @param only.nonzero  Logical: if = TRUE only the nonzero loadings are printed.
#' otherwise all loadings are printed.
#' @param perc Logical: should the loadings be standardised to unit \eqn{L_1}
#' norm (and printed as percentage contributions)?
#' @param digits Integer: number of decimal figures.
#' @param thresh Value below which loadings are considered zero and not
#' printed.
#' @param rtn Logical: should the formatted (text) table be returned?
#' @param namescomp A vector of names for the components. If NULL assigned as
#' "Comp j"
#' @param  ...  Additonal arguments for generic print, additional arguments will
#'  generate an error.
#' @return If rtn = TRUE, it returns a text table formatted as specified by the
#' arguments.
#' @note This is a wrapper for the main function in which the "dots" are disabled
#' so that only exact (or partial) prescribed arguments can be entered. 
#' @export
#' @method print spca
print.spca <- function(spca_obj, cols, only.nonzero = TRUE, perc = TRUE, digits = 3, thresh = 1E-03, rtn = FALSE, namescomp = NULL, ...){#
  goodarg <- as.list(environment())
  badarg <- eval(substitute(alist(...)))
  if (length(badarg) > 0){
    stop(paste0("\nUnused arguments: ", paste(names(badarg), collapse=", ")))
  }
# ## --
  if (is.spca(spca_obj) == TRUE){
    if (perc == TRUE){
      if(any(names(spca_obj) == "contributions"))
        A <- spca_obj$contributions
      else
        A <- make.cont(spca_obj)
    }
    else{
      if (is.null(spca_obj$As))
        A <- spca_obj$loadings
      else 
        A <- spca_obj$As
    }
  }## end is.spca == true
  else{
    if (is.matrix(spca_obj) | is.vector(spca_obj)){
      if(is.matrix(spca_obj)){
        if (perc == TRUE)
          A <- make.cont(smpc = spca_obj)
        else
          A = spca_obj
      }
      if (is.vector(spca_obj))
        if (is.vector(A))
          A <- as.matrix(A, ncol = 1)
    }    
    else
      stop(paste("The argument must be either a spca_obj object or a matrix, not a", 
                 class(spca_obj), "object"))
  } ## end if
  if (missing(cols))
    cols <- 1:ncol(A)
  else
    if (length(cols) == 1L)
      cols <- 1:cols
  if (only.nonzero == FALSE)
    rows <- 1:nrow(A)
  else{
    rows <- which(rowSums(abs(A) > thresh) > 0)
}  
  A <- as.matrix(A[rows,cols])
  ## assigns names to laodings
  if (!is.null(namescomp) & length(namescomp) == ncol(A)){
    colnames(A) <- namescomp
  }
  else{
    if (!is.null(namescomp) & length(namescomp) != ncol(A))
      message("the length of namescomp is incorrect, automatic names assigned")
    colnames(A) <- paste("Comp",1: ncol(A), sep = "")
  }
  
# # -----  formatting -
  
  if (perc == TRUE)
    fx <- format(round(A*100, max(digits-2,0)),drop0trailing = TRUE, justify = "centre")
  else
    fx <- format(round(A, digits),drop0trailing = TRUE, justify = "centre")
  names(fx) <- NULL
  nc <- nchar(fx[1L], type = "c")
  fx[abs(A) < thresh] <- paste(rep(" ", nc), collapse = "")
  #  ind = (abs(A)> thresh & abs(A) < thresh)
  #  fx[ind] = "--"
  fx <- format(fx, justify = "right" )
  if (any(class(spca_obj) == "spca")){
    vexp <- cumsum(spca_obj$vexp[cols])
    doo <- rep("-----", ifelse(is.null(ncol(fx)), 1, ncol(fx)))
    fx <- rbind(fx, doo, round(100*vexp,1))
    rownames(fx)[nrow(fx)-1] <- ""    
    rownames(fx)[nrow(fx)] <- "CVEXP"
  }  
  if (perc == TRUE)
    message("Percentage Contributions")
  else
    message("Loadings")
  if (ncol(A) == 1L){
    print(t(fx), quote = FALSE)#, ...)
  }
  else{
    print(fx, quote = FALSE)#, ...)
    cat(paste(" "))
  }
  if(rtn == TRUE){   
    return(fx)    
  }  
  else 
    invisible()
}

# #' @export
# summary <- function(spca_obj, ...){
#   UseMethod("summary", spca_obj)
# }

# #' \strong{ Prints summaries from an spca object}
# #' 
# #' Prints summaries and comparisons with the full PCA solutions for a set of LS
# #' SPCA loadings.
# #' 
# #' The summaries are printed as formatted text, if rtn = TRUE, the value
# #' returned is a numerical matrix.
# #' 
# #' For each component the following summaries are computed: \tabular{ll}{ 
# #' VEXP\tab The percentage variance explained\cr
# #' CVEXP \tab The percentage cumulative variance explained\cr
# #' RCVEXP \tab The percentage cumulative variance explained relative to that of the corresponding principal components\cr
# #' Card \tab The cardinality, that is the number of non zero loadings\cr
# #' MinLoad \tab Minimum absolute value of the non-zero loadings or contribution if perc = TRUE. }
# #' , the last row gives the minimum absolute percentage contribution, MinPContr.
# #'
# #' 
# #' @param object An spca object.
# #' @param cols A vector indicating which components should be included. Default
# #' all.  If an iteger is passed, it is set to 1:cols.
# #' @param perc Logical: should the loadings be standardised to unit L1 norm
# #' (and printed as percentage contributions)
# #' @param rtn Logical: should the summary matrix of summaries be returneded?
# #' @param prn Logical: should anything be printed? Takes priority on prnload.
# #' @param thrsehcard Value below which loadings are considered zero and not
# #' counted in the cardinality
# #' @return If rtn = TRUE, a numerical matrix with the summaries.
# #' @seealso Examples in \code{\link{SPCA}, \link{spcabe}}
# #' @export
# #' @method summary spca
# summary.spca <- function(spca_obj, cols, perc = TRUE, rtn = FALSE, prn = TRUE, 
#                         thrsehcard = 0.001){  
#                         # goodarg = as.list(environment())
#                         # badarg = eval(substitute(alist(...)))
#                         # if (length(badarg) > 0){
#                         #   stop(paste0("\nUnused arguments: ", paste(names(badarg), collapse=", ")))
#                         # }
#                         # -- 
#                         ## generic S3 method creates summaries from spca object printing for contributions
#                         if ((is.spca(spca_obj) == FALSE)){
#                           stop("summary.spca works only for spca objects")
#                         }  
#                         if( is.vector(spca_obj$loadings) )
#                           spca_obj$loadings <- as.matrix(spca_obj$loadings)
#                         if (missing(cols)){
#                           ## questo cambialo in = vexp
#                           cols <- 1:min(ncol(spca_obj$loadings), length(spca_obj$vexpPC))  
#                         }
#                         else
#                           if (length(cols) == 1L)
#                             cols <- 1:cols
#                           
#                           out <- rbind(round(spca_obj$vexp*100,1), 
#                                       round(cumsum(spca_obj$vexp)*100,1),
#                                       round(cumsum(spca_obj$vexp)/ cumsum(spca_obj$vexpPC)*100,1),
#                                       round(apply(abs(spca_obj$loadings)> thrsehcard,2,sum))
#                           )
#                           out <- as.matrix(out[,cols])
#                           colnames(out) <- paste("Comp", cols, sep = "")  
#                           rownames(out) <- 
#                             c("VEXP", "CVEXP", "RCVEXP",
#                               "Card")
#                           
#                           nc <- 4
#                           if (!is.null(spca_obj$conv)){
#                             out <- rbind(out, spca_obj$conv[cols])
#                             rownames(out)[nc + 1] <- "Converged"
#                             nc <- nc + 1
#                           }
#                           out <- rbind(out, get.minload(spca_obj, perc = perc)[cols])
#                           if (perc == FALSE)
#                             rownames(out)[nc + 1] = "MinLoad"
#                           else
#                             rownames(out)[nc + 1] = "MinCont"
#                           
#                           nc <- nc + 1            
#                           if(prn == TRUE){
#                             fx <- format(out, digits = 1, drop0trailing = FALSE, justify = "right")
#                             if (ncol(out) > 1L){ 
#                               fx[c(1:3, 5:nc),] <- apply(out[c(1:3, 5:nc),], 2, paste , "%", sep = "")  
#                               fx[4,] <- format(round(out[4,]), drop0trailing = TRUE, justify = "right",trim = TRUE, nsmall=0)  
#                               if (!is.null(spca_obj$conv)){
#                                 fx[nc -1,] <- format(round(out[nc - 1,]), drop0trailing = TRUE, justify = "right",trim = TRUE, nsmall=0)  
#                               }
#                             }
#                             else{
#                               fx[c(1:3, 5:nc),1] <- paste(out[c(1:3, 5:nc),1], "%", sep = "")  
#                               fx[4,1] <- format(round(out[4,1]), drop0trailing = TRUE, justify = "right",trim = TRUE, nsmall=0)  
#                               if (!is.null(spca_obj$conv)){
#                                 fx[nc -1,] <- format(round(out[nc - 1,]), drop0trailing = TRUE, justify = "right",trim = TRUE, nsmall=0)  
#                               }
#                               
#                             }
#                             if (perc == FALSE)
#                               fx[nc,] <- format(round(out[nc,], digits = 3), digits = 3, drop0trailing = FALSE, justify = "right")
#                             else
#                               fx[nc,] <- paste(round(out[nc,] *100,1), "%", sep = "")
# 
#                             print(fx, quote = FALSE, justify = "right")#, ...)
#                           }
#                           
#                           
#                           if (rtn == TRUE )
#                             return(out)
#                           else
#                             invisible()
# }


## plot.spca=====================
#' \strong{Plot loadings or contributions from an `spca` object}
#'
#' Plots the sparse loadings (or the corresponding percent contributions) stored in an `spca` object. Plots can be produced as linear bar plots, circular bar plots, or a
#' tile heat map. For large problems, it is recommended to set `onlynonzero = TRUE` and
#' `varnames = FALSE`.
#'
#' If `pcloadings` is provided, the plot overlays SPCA and PCA values for comparison
#' (circular bar plots with PCA are not implemented; standard bar plots are used instead).
#'
#' @param spca_obj An object of class `spca`.
#' @param nplot Integer. Number of components to plot. Defaults to `spca_obj$ncomps`.
#'   If a single integer is passed, components `1:nplot` are plotted.
#' @param plotcontributions Logical. If `TRUE`, plots contributions (typically scaled to
#'   sum to 100\% within each component); otherwise plots raw loadings.
#' @param onlynonzero Logical. If `TRUE`, plots only nonzero entries.
#' @param varnames Logical or character vector. If `TRUE`, uses row names of the loading matrix (or `VAR1, ..., VARp`. if missing). If a character vector of length `p` is supplied, it is used as the variable labels. If `FALSE`, variable labels are omitted.
#' @param vargroups Optional factor/character vector of length `p` defining groups of variables. If provided, bars/tiles are colored by group instead of by component.
#' @param plottitle Optional character. Plot title (added with `labs(title = ...)`).
#' @param stripnames Optional character vector of facet strip labels for components.
#'   Defaults to `"Comp 1"`, `"Comp 2"`, \dots, `"Comp nplot"`.
#' @param adjustLabelsCirc Optional numeric vector of length `nplot`. Additive adjustment
#'   (in degrees) to the component label angles in circular bar plots.
#' @param plottype Character. Plot type: `"bars"`, `"circular"`, or `"tiles"`. Partial
#'   matching is supported via the first letter (e.g., `"b"`, `"c"`, `"t"`).
#' @param plotgrid Logical or character. If `FALSE`, removes the background grid. If `"h"`,
#'   keeps only horizontal grid lines. If `TRUE`, uses the default grid.
#' @param legendPosition Character or numeric. Legend position. Can be one of
#'   `"bottom"`, `"right"`, `"top"`, `"left"`, or a numeric vector of length 2 giving
#'   coordinates inside the panel. For circular plots, the legend is placed on the right.
#' @param legendTitle Optional character. Legend title for the fill aesthetic.
#' @param vert Logical. If `TRUE`, flips axes for tile plots (`coord_flip()`).
#' @param pcloadings Optional numeric matrix of PCA loadings (or PCA contributions) with
#'   the same dimensions as `spca_obj$loadings`. If supplied, SPCA and PCA values are
#'   plotted together for comparison (not available for `plottype = "circular"`).
#' @param colourscale Character or character vector. One of `"cbb"`, `"printsafe"`, `"bw"`,`"ggplot"`, or a vector of hex color codes. If `vargroups` is provided, the palette is applied to groups; otherwise it is applied to the components.
#' @param returnplot Logical. If `TRUE`, returns the `ggplot2` object.
#' @param produceplot Logical. If `TRUE`, prints the plot. Useful to set `FALSE` when
#'   only returning the `ggplot2` object.
#'
#' @details
#' \code{colourscale} options: '"cbb"' is colorblind-friendly (with black), `"printsafe"` is colorblind- and printer-friendly, `"bw"` uses gray tones, and "ggplot"` uses the default ggplot2 scale. A custom vector of hex RGB colors can also be supplied. For more than 7 colours use "ggplot"
#'
#' @return If `returnplot = TRUE`, the `ggplot2` object; otherwise `NULL` (invisibly).
#' @seealso \code{\link{lsspca}}.
#' @references
#' Circular bar plot layouts follow examples from \url{https://www.r-graph-gallery.com/all-graphs/}.
#' The `cbb` palette is adapted from \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}.
#' The `printsafe` palette corresponds to `OrRd` from \url{http://colorbrewer2.org/}.
#' @export
#' @method plot spca
plot.spca <- function(spca_obj, nplot, plotcontributions = TRUE, onlynonzero = TRUE, varnames = TRUE,  vargroups = NULL, plottitle = NULL,
                      stripnames = NULL, adjustLabelsCirc = NULL,
                      plottype = c("bars", "circular", "tiles"), 
                      plotgrid = c(TRUE, "h"), legendPosition = c("bottom", "right", "top", "left"), legendTitle = NULL, vert = FALSE, pcloadings = NULL, colourscale = c("ggplot", "cbb", "printsafe", "bw"), returnplot = FALSE, produceplot = TRUE){
  
#validation=============
  p = nrow(spca_obj$laodings)
  
  if(!(any(class(spca_obj) == "spca")))
    stop("plot.spca requires an spca object as first argument")
  
  if(!(is.factor(vargroups)) && (is.vector(vargroups)))
    vargroups = vec2fac(vargroups)
  if(!is.NULL(vargroups)){
   if((!is.factor(vargroups))){
    warning("vargroups must be a character vector or a factor. Ignored.")
    vargroups = NULL
   }
    if(length(vargroups) != p){
      warning("vargroups must have length equal to the number of variables. Ignored.")
      vargroups = NULL
    }
  }
  
  if (is.character(legendPosition[1]))
    legpo = legendPosition[1]
  else
    legpo = legendPosition[1:2]
  if(grepl( "^c", plottype[1]) == TRUE){
    if (is.character(legendPosition[1]))
      legpo = "right"
    else
      legpo = c(0.5, 0.5)
  }
  thisTheme = ggplot2::theme_light() + 
    ggplot2::theme(legend.position = legpo, 
                   legend.title = element_blank(),
                   legend.text = element_text(colour="black", size = 12),
                   legend.key.size = ggplot2::unit(0.35, "cm"), 
                   legend.key.width = ggplot2::unit(0.25, "cm"),
                   panel.border = element_rect(colour = "black", linewidth = 1), 
                   panel.background = element_rect(fill = "white", 
                                                   colour = "black", linetype = 1),
                   strip.text = element_text(size = 18, color = "black"),
                   strip.background = element_rect(fill = "white", 
                                                   color = "black", linetype = 1, linewidth = ggplot2::rel(1))
    ) 
  if (plotgrid[1] == FALSE)
    thisTheme = thisTheme + 
    ggplot2::theme(panel.grid  = element_line(colour = ifelse(plotgrid[1], "grey90", NA), linewidth = 0.1))
  else 
    if (plotgrid[1] == "h")
      thisTheme = thisTheme + 
    ggplot2::theme( # remove the vertical grid lines
      panel.grid.major.x = element_blank() ,
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line( linewidth = 0.1, color = "grey90")) 
  
  if (is.spca(spca_obj) == FALSE){
    stop("plot.spca works only for spca objects")
  } 
  if(missing(nplot))
    nplot = spca_obj$ncomps
  
  if(grepl("^b", plottype[1]) == TRUE)
    plottype = "bars"
  if(grepl( "^c", plottype[1]) == TRUE)
    plottype = "circular"
  if(grepl( "^t", plottype[1]) == TRUE)
    plottype = "tiles"
  
  if((plottype == "circular") & !(is.null(pcloadings)))
    warning("\nCircular barplots with PCloadings are not implemented\n 
            (and probably too messy to be useful),\n
            using standard barplots ")
  
  
  if(grepl("^c", colourscale[1]) == TRUE)
    colourscale = "cbb"
  else
    if(grepl( "^p", colourscale[1]) == TRUE)
      colourscale = "printsafe"
  else 
    if(grepl( "^b", colourscale[1]) == TRUE)
      colourscale = "bw"
  else 
    colourscale = "ggplot"
  
  
  n <- nrow(spca_obj$loadings)
  
  plotlab = TRUE
  if(varnames[1] == TRUE){
    if (is.null(rownames(spca_obj$loadings))){
      lbl <- paste0("VAR", 1:n)
    }
    else
      lbl <- rownames(spca_obj$loadings)
  }  
  else
    if(length(varnames) == n)
      lbl <- varnames
  else{
    lbl <- paste0("VAR", 1:n)
    plotlab = NULL}
  
  
  if (is.null(stripnames))
    stripnames = paste("Comp", 1:nplot)
  else{
    if (length(stripnames) < nplot){
      warning("length of stripname must be equa to the number of plots")
      paste(stripnames, 1:nplot)
    } 
  }
  ### plots only SPCA loadings = 
  if (is.null(pcloadings)){  
    if (plotcontributions == FALSE){
      data_df <- data.frame(
        variable = factor(rep(1:n, nplot), labels = lbl),
        component = factor(rep(1:nplot, each = n), labels = stripnames),
        value = c(spca_obj$loadings[, 1:nplot])
      )
    }
    else{
      if (is.null(spca_obj$contributions))
        contributions <- scaleColsC(spca_obj$loadings[, 1:nplot], 1, rep(1, nplot))
      else 
        contributions <- spca_obj$contributions[, 1:nplot]
      data_df <- data.frame(
        variable = factor(rep(1:n, nplot), labels = lbl),
        component = factor(rep(1:nplot, each = n), labels = stripnames),
        value = c(contributions[, 1:nplot]))
    }
    ## adds factor for groups
    if (!is.null(vargroups)){
      data_df$vargroups = rep(vargroups, nplot)
    }
    
    ###  ONLY NON ZERO
    if(onlynonzero == TRUE){
      ii = (data_df$value != 0)
      data_df <- data_df[ii, ]
      # ii = (data_df$value == 0)
      # data_df[ii, ] = NA
    }
    # Create plot
    if(plottype == "circular"){
      # Set a number of 'empty bar' to add at the end of each component
      empty_bar <- 4
      to_add = data.frame( matrix(NA, empty_bar*nlevels(data_df$component), ncol(data_df)) )
      colnames(to_add) <- colnames(data_df)
      to_add$component <- rep(levels(data_df$component), each = empty_bar)
      data_df <- rbind(data_df, to_add)
      data_df <- dplyr::arrange(data_df, component)
      data_df$id <- seq(1, nrow(data_df))
      
      # Get the name and the y position of each label
      label_data <- data_df
      number_of_bar <- nrow(label_data)
      angle <- 90 - 360 * (label_data$id - 0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
      label_data$hjust <- ifelse( angle < -90, 1, 0)
      label_data$angle <- ifelse(angle < -90, angle + 180, angle)
      
      mval <- max(na.omit(data_df$value))
      label_data$pos <- data_df$value
      label_data$pos[data_df$value == 0] <- mval * (0.66)
      label_data$pos[data_df$value < 0] <- 0.05
      
      # Make the plot
      if (!is.null(vargroups)){
        pl <- ggplot2::ggplot(data_df, aes(x = as.factor(id), y = value, fill = vargroups))
      }
      else{
        pl <- ggplot2::ggplot(data_df, aes(x = as.factor(id), y = value, fill = component)) 
      }
      pl = pl + ggplot2::geom_bar(stat = "identity", alpha = 0.75, 
                                  color = ifelse(colourscale[1] == "printsafe", "black", NA)) +
        ggplot2::ylim(-0.5 - median(abs(na.omit(data_df$value))), mval + 0.2) + 
        ggplot2::theme_light() +
        ggplot2::theme(
          legend.position = legpo, #"right"
          axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.margin = ggplot2::unit(rep(-1,4), "cm") 
        ) +  labs(fill = legendTitle) +
        ggplot2::geom_abline(slope = 0, intercept = 0) 
      
      #      if(!is.null(vargroups)){-
      indg = (!is.na(data_df$variable))
      minloads = tapply(X = as.numeric(data_df$id[indg]), INDEX = data_df$component[indg], min) 
      # maxloads = tapply(X = as.numeric(data_df$id[indg]), INDEX = data_df$component[indg], max) 
      medloads = tapply(data_df$id[indg],  INDEX = data_df$component[indg], median)
      # medangle = tapply(angle[indg],  INDEX = data_df$component[indg], median)
      lia = min(data_df$value[indg]) - 0.05
      
      medangle = 360*medloads/number_of_bar
      
      # anglela <- 90 - 360  * (medloads - 0.5)/nplot #+  
      #   c(-90, -270, -270, 270)[1:nplot]
      #anglela <- medangle #+ ifelse(((medangle < -100) & (medangle > -270)), 180, 90)
      anglela = 180 - medangle #- (1:nplot)* 45
      anglela[((anglela >90) | (anglela < -90))] <- anglela[((anglela >90) | (anglela < -90))] + 180
      if(!is.null(adjustLabelsCirc)){
        if(length(adjustLabelsCirc) == nplot)
          anglela = anglela + adjustLabelsCirc
        else
          warning("need to pass as many values to adjustLabelCirc as nplot")
      }      
      
      pl = pl + ggplot2::annotate("text", x = medloads, y = rep(lia, nlevels(data_df$component)), 
                                  label = paste("comp", 1:nplot), hjust = 0.5, size = 4, angle = anglela,
                                  fontface = "bold") + ggplot2::coord_polar()
      #}
      
      if (!is.null(plotlab)) 
        pl <- pl + ggplot2::geom_text(data = label_data, aes(x = id, y = pos + 0.05, label = variable, hjust = hjust), 
                                      color = "black", fontface = "bold", alpha = 0.85, size = 3.5, 
                                      angle = label_data$angle, inherit.aes  =  FALSE ) 
    }
##plot bars=======================    
    else{ 
      if(plottype == "bars"){
        # linear barplots
        nrows <- ceiling(nplot/3)
        ncols <- ceiling(nplot/nrows)
        if (!is.null(vargroups))
          pl <- ggplot2::ggplot(data_df, aes(x = variable, y = value, fill = vargroups)) 
        else        
          pl <- ggplot2::ggplot(data_df, aes(x = variable, y = value, fill = component)) 
       pl = pl + ggplot2::geom_bar(stat = "Identity", color = ifelse(colourscale[1] == "printsafe", "black", NA)) +
          ggplot2::facet_wrap(facets = vars(component), ncol = ncols, nrow = nrows) +
          ggplot2::geom_abline(slope = 0, intercept = 0) + thisTheme +
          ggplot2::xlab("variables") + ggplot2::ylab( ifelse(plotcontributions == TRUE, "contributions", "loadings"))
        if (!is.null(plotlab))
          pl <- pl + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)) 
        else
          pl <- pl +  ggplot2::theme(axis.ticks = element_blank(), axis.text.x = element_blank())
       if(plotcontributions)
        pl =  pl + ggplot2::scale_y_continuous(labels = scales::percent)
       
        if (!is.null(vargroups))
          pl <- pl +  ggplot2::theme(legend.position = legpo)#"right")
      }
      else
        if (plottype == "tiles"){
          ## this is colorbrewer red2Blue_pal
          tile_pal = c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", 
                       "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
          
          pl = ggplot2::ggplot(data_df, aes(variable, component)) +
            geom_tile(aes(fill = value), colour = "gray75") +
            theme_bw() + scale_fill_gradientn(colours = tile_pal, limits =c(-1, 1),
                                              name = ifelse(plotcontributions, "Contributions", "Loadings")) + 
            theme(legend.position = legpo)#"bottom") 
          pl <- pl  + ggplot2::geom_abline(intercept = (1:nplot) + 0.5, slope = 0, colour = "grey75") +
            ggplot2::geom_vline(xintercept = (1:length(unique(data_df$variable))) + 0.5, colour = "grey75") + 
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                  axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8))
          
          if (vert == TRUE)
            pl <- pl + coord_flip()
          
        }
    }
  }
  # plot with Pc's loadings =============#
  else{ 
    if (plottype != "tiles"){
      cs <- apply(pcloadings< 0, 2, all)
      pcloadings[, cs] <- -pcloadings[, cs]
      
      for (i in 1:nplot){
        ii <- spca_obj$ind[[i]][1]
        si <- sign(spca_obj$loadings[ii, i] * pcloadings[ii, i])
        if (si == -1){
          message(paste0("changed sign to loadings component ", i ))
          spca_obj$loadings[, i]  <- -spca_obj$loadings[, i]
          if (!is.null(spca_obj$contributions)) 
            spca_obj$loadlist[[i]] <- -spca_obj$loadlist[[i]] 
        }
      }   
      if (plotcontributions == FALSE){
        # if(onlynonzero == TRUE){
        #   ii = (pcloadings != 0)
        #   pcloadings <- pcloadings[ii, , drop = FALSE]
        # }
        
        data_df <- data.frame(
          variable = factor(c(unlist(spca_obj$ind[1:nplot]), rep(1:n, nplot)), labels = lbl),
          component = factor(c(rep(1:nplot, times = spca_obj$card[1:nplot]), rep(1:nplot, each = n)), 
                             labels = stripnames),
          value = c(unlist(spca_obj$loadlist[1:nplot]), c(pcloadings[, 1:nplot]))
        )
      }
      else{
        contributions <- unlist(
          sapply(spca_obj$loadlist[1:nplot], function(x) x/sum(abs(x)))
        )
        pccontributions <- scaleColsC(pcloadings[, 1:nplot], 1, rep(1, nplot))
        # if(onlynonzero == TRUE){
        #   ii = (pccontributions != 0)
        #   pccontributions <- pccontributions[ii, , drop = FALSE]
        # }
        
        data_df = data.frame(
          variable = factor(c(unlist(spca_obj$ind[1:nplot]), rep(1:n, nplot)), labels = lbl),
          component = factor(c(rep(1:nplot, times = spca_obj$card[1:nplot]), rep(1:nplot, each = n)), 
                             labels = stripnames),
          value = c(contributions, c(pccontributions[, 1:nplot]))
        )
      }
      data_df$method <- factor(c(rep(1, sum(spca_obj$card[1:nplot])),
                                 rep(2, n * nplot)), labels = c("SPCA", "PCA"))
      nrows <- ceiling(nplot/3)
      ncols <- ceiling(nplot/nrows)
      
      pl <- ggplot2::ggplot(data_df, aes(x = variable, y = value, fill = method)) + 
        ggplot2::geom_bar(stat = "Identity", position = position_dodge(), 
                          color = ifelse(colourscale[1] == "printsafe", "black", NA)) +
        ggplot2::facet_wrap(facets = vars(component), ncol = ncols, nrow = nrows) + 
        ggplot2::theme_classic() +
        ggplot2::theme(
          legend.position = legpo, #"right",
          legend.title=element_blank(),
          # axis.text = element_blank(),
          # axis.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black", linetype = 1)
          #,plot.margin = ggplot2::unit(rep(-1,4), "cm") 
        ) +
        ggplot2::xlab("variables") + 
        ggplot2::ylab( ifelse(plotcontributions == TRUE, "contributions", "loadings")) 
      if (!is.null(plotlab))
        pl <- pl + ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)) 
      else
        pl <- pl + ggplot2::theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    }
    else{ ## plot tiles with pcs
      
      rownames(pcloadings) = lbl 
      if(plotcontributions == TRUE){
        value = c(scaleColsC(spca_obj$loadings[, 1:nplot], 1, rep(1, nplot)))
        PCvalue = c(scaleColsC(pcloadings[, 1:nplot], 1, rep(1, nplot)))
      }
      else{
        value = c(spca_obj$loadings[, 1:nplot])
        PCPCvalue = c(pcloadings[, 1:nplot]) 
      }
      SPC_df = data.frame(Variable = factor(rep(1:length(lbl), nplot), labels = lbl),
                          Component = factor(rep(1:nplot, each = length(lbl)), 
                                             labels = stripnames),
                          value = value,
                          varNum = rep(1:length(lbl), nplot),
                          compNum = rep(1:nplot, each = length(lbl)))
      
      PC_df = data.frame(Variable = factor(rep(1:length(lbl), nplot), labels = lbl),
                         Component = factor(rep(1:nplot, each = length(lbl)), 
                                            labels = stripnames),
                         value = PCvalue,
                         varNum = rep(1:length(lbl), nplot),
                         compNum = rep(1:nplot, each = length(lbl)) + 0.5)
      
      data_df = rbind(SPC_df, PC_df)
      
      lab_y = factor(1:(2 * nplot), labels = paste(c("Comp", "PC"), rep(1:nplot, each = 2)))
      
      tile_pal = c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#F7F7F7", 
                   "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
      
      pl = ggplot(data_df, aes(xmin = varNum, xmax = varNum + 1, 
                               ymin = compNum, ymax = compNum + 0.5, fill = value)) +
        geom_rect() +   theme_bw() +
        scale_fill_gradientn(colours = tile_pal, limits =c(-1, 1)) +
        scale_x_continuous(breaks = seq(1.5, length(lbl) + 0.5, 1), 
                           labels = lbl, expand = c(0,0)) +
        scale_y_continuous(breaks = seq(1.25, (length(lab_y)/2) + 1, 0.5), 
                           labels = lab_y, expand = c(0,0)) + 
        geom_abline(intercept = (1:nplot) + 1, slope = 0, colour = "black", size = 1.5) +
        geom_abline(intercept = (1:nplot) + 0.5, slope = 0, colour = "gray75", size = 1) +
        geom_vline(xintercept = (1:length(lbl)) + 1, colour = "grey75") +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.ontop = TRUE,
              panel.background = element_rect(fill = "transparent"),
              axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8)) 
    }
  }
  if (plottype != "tiles"){
    if (colourscale[1] == "cbb"){{## use colourblind friendly with black
      cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      if (!is.null(pcloadings))
        cbb <- cbb[c(1,7)]
      
      if(!is.null(vargroups)){
        pl <- pl  + ggplot2::scale_fill_manual(limits = levels(data_df$vargroups), values = cbb) 
      }
      
      else
        pl <- pl  + ggplot2::scale_fill_manual(values = cbb) 
      colourscale <- "ggplot"
    }
    }  
    
    if (colourscale[1] == "printsafe"){## use RColorBrewer::brewer.pal("Colors in OrRd", "OrRd")
      pf.pal <- 
        rev(c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
      if (!is.null(pcloadings))
        pf.pal = pf.pal[c(1, 4)]
      if(!is.null(vargroups)){
        pl <- pl  + ggplot2::scale_fill_manual(limits = levels(data_df$vargroups), values = pf.pal) 
      }
      # }
      else 
        pl <- pl  + ggplot2::scale_fill_manual(values = pf.pal) 
      colourscale <- "ggplot"
    }
    if ((colourscale[1] == "bw")){
      pl <- pl  + ggplot2::scale_fill_grey() 
      colourscale <- "ggplot"
    }
    else{
      if(colourscale[1] != "ggplot"){
        if(!is.null(vargroups)){
          pl <- pl  + ggplot2::scale_fill_manual(limits = levels(data_df$vargroups),
                                                 values = colourscale[1:nplot][ii])
        }
        else{
          pl <- pl  + ggplot2::scale_fill_manual(values = colourscale[1:nplot][ii])
        }
      }
    }
  }
  if (!is.null(plottitle))
    pl = pl + labs(title = plottitle)
  if(produceplot == TRUE)
    print(pl)
  if(returnplot == TRUE)
    return(pl)
  else
    invisible()
}

## showload ==================
#' \strong{Shows the non-zero loadings separately for each component.}
#' 
#' Useful for large matrices to see the loadings at the same 
#' time or to assign long descriptive names.
#' 
#' @param spca_obj A list of spca objects, typically from SPCA.  It
#' can also be a simple matrix of loadings.
#' @param cols A vector containg the indices of the loadings to be shown.  Can be a single value. if missing all loadings are shown: If an integer is
#' passed, only that dimension will be returned.
#' @param perc Logical: should the loodings be standardised to unit \eqn{L_1} norm (and printed as percentage contributions).
#' @param digits Number of decimal digits to show.
#' @param  variablesnames Hybrid: if not FALSE, need to pass a vector of varaiable names.
#' @param thresh Loadings with absolute value below this are considered zero.
#' @param rtn Logical: should the text table of loadings and the matrix of summaries be returneded?
#' @details variablesnames must have the names of the p variables in the first p
#' positions.  
#' @return If rtn = TRUE, it returns a list with the loadings.
#' @seealso \link{PMlsspca-package}.
#' @export 
showload = function(spca_obj, cols, perc = TRUE, digits = 3,  variablesnames = FALSE, 
                    thresh = 0.001, rtn = FALSE){
  
## function that prints nonzero loadings one component at the time from an spca object
  
  if(missing(cols)){
    cols = 1:ncol(spca_obj$loadings)
  }
  if (!any(class(spca_obj) == "spca")){
    if (is.matrix(spca_obj) | is.vector(spca_obj))
      A = as.matrix(spca_obj)[,cols]
    else
      stop("need an spca object or an array of loadings")
  }
  else{
    A = as.matrix(spca_obj$loadings[,cols])
    if (length(cols) == 1)
      rownames(A) =  rownames(spca_obj$loadings)
  }
  if (perc == TRUE){
    dimna = dimnames(A)
    A = scaleColsC(A, 1, rep(1, ncol(A)))
    dimnames(A) = dimna
  }
  loads = list()
  if (! is.null(rownames(A)))
    variablesnames = rownames(A)
  if (is.factor(variablesnames))
    variablesnames = as.character(variablesnames)
  if (is.vector( variablesnames) ){ 
    if (length( variablesnames) < nrow(A))
      stop(" variablesnames must be a vector of length equal to number of loadings")
    else
      rownames(A) =  variablesnames[1:nrow(A)]
  }
  if (perc == TRUE)
    message("Percent Contributions")
  else
    message("Norm 1 Loadings")
  for(i in 1:length(cols)){
    loads[[i]] = A[abs(A[,i])> thresh,i]
    print(paste("Component", cols[i]))
    if (perc == TRUE){
      a = paste(round(100 * loads[[i]], max(0,digits-2)), "%", sep = "")
      names(a) = names(loads[[i]])
      print(a, quote = FALSE, justify = "left")
      writeLines(" ")#paste(, quote = FALSE)
    }
    else{
      print(loads[[i]], digits = digits, justify = "left")
      writeLines(" ")
    }
  }
  if (rtn == TRUE)
    return(loadings)
  else
    invisible()
  
}

#new.spca==================
#' \strong{Constructs an object of class spca from a set of loadings}
#' 
#' @param A real matrix, a matrix of loadings.
#' @param S the variance or correlation matrix from which the loadings where computed.
#' @param X real matrix, the data matrix from which the loadings where computed, optional
#' @param method string, the name of the method used to compute the loadings, optional.
#' @return An object is of class spca with several added components, 
#' FALSE otherwise.
#' @details The \code{spca} object can be used as argument to the spca methods.
#'  It contains \describe{ 
#'  \item{loadings}{the matrix of loadings.}
#'  \item{contributions}{the matrix of loadings scaled to unit L1 norm.} 
#'  \item{ncomps}{the number of components.}
#'  \item{cardinality}{the vector with the number of nonzero loadings in 
#'  each each component.}
#'  \item{ind}{a list with the indices of the variables corresponding the nonzero loadings.}
#'  \item{vexp}{a vector of percent variance explained by the components.} 
#'  \item{cvexp}{a vector of percent cumulative variance explained by the components.}
#'  \item{vexpPC}{a vector of variance explained by the PCs.}
#'  \item{rvexp}{a vector of proportion of cumulative variance explained over that explained by the 
#'  corresponding PCs.}
#'  \item{loadlist}{a list with only the nonzero loadings.}
#'  \item{scores}{the components' scores, if X is passed.} 
#'  \item{corComp}{A matrix of the correlations between components.}
#'  \item{vif}{a list of variance inflation factors for the variables in each components, 
#'  if X is passed.}
#'  \item{method}{the string passed as such, if method is passed.}
#' }
#' @seealso \code{\link{PMlsspca-package}} for a description of an \code{spca} object.
#' @export 
new_spca = function(A, S, X = NULL, method = NULL){
  
#            validation   ==============P
  if(is.data.frame(A)) A = as.matrix(A)
  if (!is.matrix(A)){
    stop("A must be a matrix of loadings")
  }
  if((missing(S) && !(is.null(X))))
  if (!(is.matrix(A) & isSymmetric(S))){
    stop("S must be a variance or correlation matrix")
  }
  
  obj = list()
  obj$loadings = A
  obj$contributions = scaleColsC(A, 1, rep(1, ncol(A)))
  dimnames(obj$contributions) = dimnames(A)
  obj$ncomps = ncol(A)
  obj$cardinality = colSums(A != 0)
  
  obj$ind = list()
  loadlist = list()
  for(i in 1:obj$ncomps){
    obj$ind[[i]] = which(A[, i] != 0)
    loadlist[[i]] = A[obj$ind[[i]], i]
    
    names(obj$ind[[i]]) = rownames(A)[obj$ind[[i]]]
    names(loadlist[[i]]) = rownames(A)[obj$ind[[i]]]
  }
  
  
  aa = makeVexpSC(obj$loadings, S)  
  s_ee = eigen(S)
  totv = sum(s_ee$values)
  obj$vexp = aa$vexp/totv
  obj$cvexp = aa$cvexp/totv
  
  obj$vexpPC = s_ee$values[1:obj$ncomps]/totv
  obj$rcvexp = obj$cvexp/cumsum(obj$vexpPC)
  obj$loadlist = loadlist
  
  if(!is.null(X)){
    obj$scores = X %*% A
    obj$corComp = cor(obj$scores)
  }
  else
    obj$corComp = MkCorCompMat(A, S)
  class(obj) = c("list", "spca")
  obj$vif <- make_vif.spca(obj, S, prn = FALSE)
  if(!is.null(method))  
    obj$method = method
  return(obj)
}

#aggregate_by_scale========
#' \strong{Aggregate loadings or contributions by scale}
#'
#' Computes scale-level sums by aggregating an input vector or the columns of an
#' input matrix/data.frame according to a grouping index (e.g., scale membership
#' of variables).
#'
#' @param x A numeric vector, or a numeric matrix/data.frame containing loadings or contributions.
#' @param ind A vector or factor of length equal to `length(x)` (if `x` is a vector) or `nrow(x)` (if `x` is a matrix/data.frame), giving the scale/group label for each variable.
#' @param only.nonzero Logical. If `TRUE`, scales with zero total (based on
#'   `rowSums(abs(out))`) are removed from the output.
#'
#' @return A numeric vector (if `x` is a vector) or numeric matrix (if `x` is a
#'   matrix/data.frame) of sums aggregated by scale. Rows correspond to scales
#'   and columns correspond to columns of `x`.
#'
#' @export
aggregate_by_scale <- function(x, ind, only.nonzero = TRUE) {
  if (is.vector(x)) {
    out = tapply(x, ind, sum)
  } else {
    out = apply(x, 2, function(y, ii) tapply(y, ii, sum), ii = ind)
    if (only.nonzero)
      out = droplevels(out[rowSums(abs(out)) > 0, ])
  }
  return(out)
}

#summary.spca==============
#' \strong{ Prints summaries from an spca object}
#' 
#' Prints summaries and comparisons with the full PCA solutions for a set of LS
#' SPCA loadings.
#' 
#' The summaries are printed as formatted text, if rtn = TRUE, the value
#' returned is a numerical matrix.
#' 
#' For each component the following summaries are computed: \tabular{ll}{ 
#' VEXP\tab The percentage variance explained\cr
#' CVEXP \tab The percentage cumulative variance explained\cr
#' RVEXP \tab The variance explained relative to the corresponding PC\cr
#' RCVEXP \tab The cumulative variance explained relative to the corresponding PCs\cr
#' Card \tab The cardinality, that is the number of non zero loadings\cr
#' MinLoad/MinCont \tab Minimum absolute value of the non-zero loadings or contribution }
#'
#' @param spca_obj An spca object.
#' @param cols A vector indicating which components should be included. Default
#' all. If an integer is passed, it is set to 1:cols.
#' @param contribution Logical: should the loadings be standardised to unit L1 norm
#' (and printed as percentage contributions)
#' @param variance_metrics Character vector: which variance metrics to include.
#' Options: "relative" (RVEXP), "cumulative_relative" (RCVEXP), "both", or "none".
#' Default is "cumulative_relative".
#' @param minload Logical: should minimum loading/contribution be included?
#' @param rtn Logical: should the summary matrix be returned?
#' @param prn Logical: should anything be printed?
#' @param thrsehcard Value below which loadings are considered zero and not
#' counted in the cardinality
#' @return If rtn = TRUE, a numerical matrix with the summaries.
#' @seealso Examples in \code{\link{lsspca}}.
#' @export
#' @method summary spca
summary.spca <- function(spca_obj, cols, contribution = TRUE, 
                         variance_metrics = "cumulative_relative",
                         minload = FALSE, rtn = FALSE, prn = TRUE, 
                         thrsehcard = 0.001) {
  
  # Validation
  if (!is.spca(spca_obj)) {
    stop("summary.spca works only for spca objects")
  }
  
  if (is.vector(spca_obj$loadings)) {
    spca_obj$loadings <- as.matrix(spca_obj$loadings)
  }
  
  # Determine columns
  if (missing(cols)) {
    cols <- 1:min(ncol(spca_obj$loadings), length(spca_obj$vexpPC))
  } else if (length(cols) == 1L) {
    cols <- 1:cols
  }
  
  # Validate variance_metrics
  valid_metrics <- c("relative", "cumulative_relative", "both", "none")
  if (!variance_metrics %in% valid_metrics) {
    stop("variance_metrics must be one of: ", paste(valid_metrics, collapse = ", "))
  }
  
  # Build summary matrix
  out <- rbind(
    VEXP = spca_obj$vexp * 100,
    CVEXP = cumsum(spca_obj$vexp) * 100
  )
  
  # Add variance comparison metrics based on user choice
  if (variance_metrics %in% c("relative", "both")) {
    out <- rbind(out, RVEXP = spca_obj$vexp / spca_obj$vexpPC * 100)
  }
  
  if (variance_metrics %in% c("cumulative_relative", "both")) {
    out <- rbind(out, RCVEXP = cumsum(spca_obj$vexp) / cumsum(spca_obj$vexpPC) * 100)
  }
  
  # Add cardinality
  out <- rbind(out, Card = apply(abs(spca_obj$loadings) > thrsehcard, 2, sum))
  
  # Add convergence if available
  if (!is.null(spca_obj$conv)) {
    out <- rbind(out, Converged = spca_obj$conv)
  }
  
  # Add minimum loading/contribution if requested
  if (minload) {
    if (contribution) {
      min_vals <- get.minload(spca_obj, perc = TRUE) * 100
      min_label <- "MinCont"
    } else {
      min_vals <- get.minload(spca_obj, perc = FALSE)
      min_label <- "MinLoad"
    }
    out <- rbind(out, temp = min_vals)
    rownames(out)[nrow(out)] <- min_label
  }
  
  # Select requested columns
  out <- as.matrix(out[, cols, drop = FALSE])
  colnames(out) <- paste0("Comp", cols)
  
  # Print formatted output
  if (prn) {
    out_formatted <- format_summary_matrix(out, contribution)
    print(out_formatted, quote = FALSE, justify = "right")
  }
  
  # Return if requested
  if (rtn) {
    return(out)
  } else {
    invisible()
  }
}


# Format summary matrix for printing
# @keywords internal
format_summary_matrix <- function(out, contribution) {
  
  # Define which rows get which formatting
  percentage_rows <- c("VEXP", "CVEXP", "RVEXP", "RCVEXP", "MinCont")
  
  integer_rows <- c("Card", "Converged")
  decimal_rows <- "MinLoad"
  
  # Format each row
  fx <- matrix("", nrow = nrow(out), ncol = ncol(out))
  rownames(fx) <- rownames(out)
  colnames(fx) <- colnames(out)
  
  for (i in 1:nrow(out)) {
    row_name <- rownames(out)[i]
    
    if (row_name %in% percentage_rows) {
      # Format as percentage with 1 decimal
      fx[i, ] <- paste0(format(round(out[i, ], 1), nsmall = 1, 
                               drop0trailing = FALSE, justify = "right"), "%")
      
    } else if (row_name %in% integer_rows) {
      # Format as integer
      fx[i, ] <- format(round(out[i, ]), drop0trailing = TRUE, 
                        justify = "right", trim = TRUE, nsmall = 0)
      
    } else if (row_name %in% decimal_rows) {
      # Format with 3 decimals (not percentage)
      fx[i, ] <- format(round(out[i, ], 3), nsmall = 3, 
                        drop0trailing = FALSE, justify = "right")
    }
  }
  
  return(fx)
}

