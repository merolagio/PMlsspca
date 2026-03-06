
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
#' @param x An spca object.
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
#' @param ... Further arguments; currently ignored.
#' @return If rtn = TRUE, it returns a text table formatted as specified by the
#' arguments.
#' @note This is a wrapper for the main function in which the "dots" are disabled
#' so that only exact (or partial) prescribed arguments can be entered. 
#' @export
#' @method print spca
print.spca <- function(x, cols, only.nonzero = TRUE, perc = TRUE, digits = 3, thresh = 1E-03, rtn = FALSE, namescomp = NULL, ...){#
  goodarg <- as.list(environment())
  badarg <- eval(substitute(alist(...)))
  if (length(badarg) > 0){
    stop(paste0("\nUnused arguments: ", paste(names(badarg), collapse=", ")))
  }
  # ## --
  if (is.spca(x) == TRUE){
    if (perc == TRUE){
      if(any(names(x) == "contributions"))
        A <- x$contributions
      else
        A <- make.cont(x)
    }
    else{
      if (is.null(x$As))
        A <- x$loadings
      else 
        A <- x$As
    }
  }## end is.spca == true
  else{
    if (is.matrix(x) | is.vector(x)){
      if(is.matrix(x)){
        if (perc == TRUE)
          A <- make.cont(smpc = x)
        else
          A = x
      }
      if (is.vector(x))
        if (is.vector(A))
          A <- as.matrix(A, ncol = 1)
    }    
    else
      stop(paste("The argument must be either a x object or a matrix, not a", 
                 class(x), "object"))
  } ## end if
  if (missing(cols) || is.null(cols))
    cols <- 1:ncol(A)
  else
    if (length(cols) == 1L)
      cols <- 1:cols
  if (only.nonzero == FALSE)
    rows <- 1:nrow(A)
  else{
    rows <- which(rowSums(abs(A[, cols]) > thresh) > 0)
  }  
  A <- as.matrix(A[rows, cols])
  ## assigns names to laodings
  if ((!is.null(namescomp)) & (length(namescomp) >= length(cols))){
    colnames(A) <- namescomp[cols]
  }
  else{
    if (!(is.null(namescomp)) & (length(namescomp) < length(cols)))
      message("the length of namescomp is incorrect, automatic names assigned")
    colnames(A) <- paste("Comp",cols, sep = "")
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
  if (any(class(x) == "spca")){
    vexp <- cumsum(x$vexp[cols])
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


# plot.spca=====================
#' \strong{Plot loadings or contributions from an `spca` object}
#'
#' Plots the sparse loadings (or the corresponding percent contributions) stored in an `spca` object. Plots can be produced as linear bar plots, circular bar plots, or a
#' tile heat map. For large problems, it is recommended to set `onlynonzero = TRUE` and
#' `varnames = FALSE`.
#'
#' If `pcloadings` is provided, the plot overlays SPCA and PCA values for comparison
#' (circular bar plots with PCA are not implemented; standard bar plots are used instead).
#'
#' @param x An object of class `spca`.
#' @param nplot Integer. Number of components to plot. Defaults to `x$ncomps`.
#'   If a single integer is passed, components `1:nplot` are plotted.
#' @param plotcontributions Logical. If `TRUE`, plots contributions (typically scaled to
#'   sum to 100\% within each component); otherwise plots raw loadings.
#' @param onlynonzero Logical. If `TRUE`, plots only nonzero entries.
#' @param varnames Logical or character vector. If `TRUE`, uses row names of the loading matrix (or `VAR1, ..., VARp`. if missing). If a character vector of length `p` is supplied, it is used as the variable labels. If `FALSE`, variable labels are omitted.
#' @param vargroups Optional factor/character vector of length `p` defining groups of variables. If provided, bars/tiles are colored by group instead of by component.
#' @param plottitle Optional character. Plot title (added with `labs(title = ...)`).
#' @param x_axis_lab label for x axis. default is "variables"
#' @param stripnames Optional character vector of facet strip labels for components.
#'   Defaults to `"Comp 1"`, `"Comp 2"`, \dots, `"Comp nplot"`.
#' @param adjustLabelsCirc Optional numeric vector of length `nplot`. Additive adjustment
#'   (in degrees) to the component label angles in circular bar plots.
#' @param plottype Character. Plot type: `"bars"`, `"circular"`, or `"tiles"`. Partial
#'   matching is supported via the first letter (e.g., `"b"`, `"c"`, `"t"`).
#' @param plotgrid Hybrid logical or character. If `FALSE`, removes the background grid. If `"h"`,
#'   keeps only horizontal grid lines. If `TRUE`, uses the default grid.
#' @param legendPosition Character or numeric. Legend position. Can be one of
#'   `"bottom"`, `"right"`, `"top"`, `"left"`, or a numeric vector of length 2 giving
#'   coordinates inside the panel. For circular plots, the legend is placed on the right.
#' @param legendTitle Optional character. Legend title for the fill aesthetic.
#' @param vert Logical. If `TRUE`, flips axes for tile plots (`coord_flip()`).
#' @param pcloadings Optional numeric matrix of PCA loadings (or PCA contributions) with the same dimensions as `x$loadings`. If supplied, SPCA and PCA values are
#'   plotted together for comparison (not available for `plottype = "circular"`).
#' @param colourscale Character or character vector. One of `"cbb"`, `"printsafe"`, `"bw"`,`"ggplot"`, or a vector of hex color codes. If `vargroups` is provided, the palette is applied to groups; otherwise it is applied to the components.
#' @param returnplot Logical. If `TRUE`, returns the `ggplot2` object.
#' @param produceplot Logical. If `TRUE`, prints the plot. Useful to set `FALSE` when only returning the `ggplot2` object.
#' @param ... Further arguments; currently ignored.
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
plot.spca <- function(x, nplot, plotcontributions = TRUE, onlynonzero = TRUE, varnames = TRUE,  vargroups = NULL, plottitle = NULL, stripnames = NULL, x_axis_lab = "variables", adjustLabelsCirc = NULL, plottype = c("bars", "circular", "tiles"),  plotgrid = c("h", TRUE, NA), legendPosition = c("bottom", "right", "top", "left"), legendTitle = NULL, vert = FALSE, pcloadings = NULL, colourscale = c("ggplot", "cbb", "printsafe", "bw"), returnplot = FALSE, produceplot = TRUE, ...){
  
#validation=============
  p = nrow(x$loadings)
  
  if(!(any(class(x) == "spca")))
    stop("plot.spca requires an spca object as first argument")
  
  if(!(is.factor(vargroups)) && (is.vector(vargroups)))
    vargroups = vec2fac(vargroups)
  if(!is.null(vargroups)){
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
  
  if (is.spca(x) == FALSE){
    stop("plot.spca works only for spca objects")
  } 
  if(missing(nplot))
    nplot = x$ncomps
  
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
  
  
  n <- nrow(x$loadings)
  
  plotlab = TRUE
  if(varnames[1] == TRUE){
    if (is.null(rownames(x$loadings))){
      lbl <- paste0("VAR", 1:n)
    }
    else
      lbl <- rownames(x$loadings)
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
        value = c(x$loadings[, 1:nplot])
      )
    }
    else{
      if (is.null(x$contributions))
        contributions <- scaleColsC(x$loadings[, 1:nplot], 1, rep(1, nplot))
      else 
        contributions <- x$contributions[, 1:nplot]
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
      data_df <- data_df[order(data_df$component), ]
      data_df$id <- seq(1, nrow(data_df))
      
      # Get the name and the y position of each label
      label_data <- data_df
      number_of_bar <- nrow(label_data)
      angle <- 90 - 360 * (label_data$id - 0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
      label_data$hjust <- ifelse( angle < -90, 1, 0)
      label_data$angle <- ifelse(angle < -90, angle + 180, angle)
      
      mval <- max(stats::na.omit(data_df$value))
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
        ggplot2::ylim(-0.5 - stats::median(abs(stats::na.omit(data_df$value))), mval + 0.2) + 
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
      medloads = tapply(data_df$id[indg],  INDEX = data_df$component[indg], stats::median)
      # medangle = tapply(angle[indg],  INDEX = data_df$component[indg], stats::median)
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
          ggplot2::xlab(x_axis_lab) + ggplot2::ylab( ifelse(plotcontributions == TRUE, "contributions", "loadings"))
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
        ii <- x$ind[[i]][1]
        si <- sign(x$loadings[ii, i] * pcloadings[ii, i])
        if (si == -1){
          message(paste0("changed sign to loadings component ", i ))
          x$loadings[, i]  <- -x$loadings[, i]
          if (!is.null(x$contributions)) 
            x$loadlist[[i]] <- -x$loadlist[[i]] 
        }
      }   
      if (plotcontributions == FALSE){
        # if(onlynonzero == TRUE){
        #   ii = (pcloadings != 0)
        #   pcloadings <- pcloadings[ii, , drop = FALSE]
        # }
        
        data_df <- data.frame(
          variable = factor(c(unlist(x$ind[1:nplot]), rep(1:n, nplot)), labels = lbl),
          component = factor(c(rep(1:nplot, times = x$card[1:nplot]), rep(1:nplot, each = n)), 
                             labels = stripnames),
          value = c(unlist(x$loadlist[1:nplot]), c(pcloadings[, 1:nplot]))
        )
      }
      else{
        contributions <- unlist(
          sapply(x$loadlist[1:nplot], function(x) x/sum(abs(x)))
        )
        pccontributions <- scaleColsC(pcloadings[, 1:nplot], 1, rep(1, nplot))
        # if(onlynonzero == TRUE){
        #   ii = (pccontributions != 0)
        #   pccontributions <- pccontributions[ii, , drop = FALSE]
        # }
        
        data_df = data.frame(
          variable = factor(c(unlist(x$ind[1:nplot]), rep(1:n, nplot)), labels = lbl),
          component = factor(c(rep(1:nplot, times = x$card[1:nplot]), rep(1:nplot, each = n)), 
                             labels = stripnames),
          value = c(contributions, c(pccontributions[, 1:nplot]))
        )
      }
      data_df$method <- factor(c(rep(1, sum(x$card[1:nplot])),
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
        value = c(scaleColsC(x$loadings[, 1:nplot], 1, rep(1, nplot)))
        PCvalue = c(scaleColsC(pcloadings[, 1:nplot], 1, rep(1, nplot)))
      }
      else{
        value = c(x$loadings[, 1:nplot])
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
  obj$vif <- make_vif.spca(obj, S)
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
#' @param perc if TRUE the aggregated scale contributions are returned as percentages.
#' @param digits If integer, the number of digits to which to round the aggregated values. If FALSE no rounding is done.  
#'
#' @details
#' If loadings are passed and perc == TRUE, perc is turned to FALSE and digits to 3 because it is not appropriate to express L2 unit loadings as percentages. 
#' 
#' @return A numeric vector (if `x` is a vector) or numeric matrix (if `x` is a
#'   matrix/data.frame) of sums aggregated by scale. Rows correspond to scales
#'   and columns correspond to columns of `x`.
#'
#' @export
aggregate_by_scale <- function(x, ind, only.nonzero = TRUE, perc = T, digits = ifelse(perc, 1, 3)) {
  #vector
if (is.vector(x)) {
  if(length(ind) != length((x)))
    stop("contributions and ind must have the same length")
  else{
    out = tapply(x, ind, sum)
    if ((abs(sum(abs(x)) - 1) > 10e-3) && perc == T){
      warning("turning the sum of loadings to percentage is not appropriate")
      perc = FALSE
      digits = 3
    }
    if(only.nonzero)
      out = out[abs(out) > 10e-5]
    }
  } 
 # matrix
  else {
    if(length(ind) != nrow((x)))
      stop("The length of ind must be the same as the rows of x ")
    else{
      if (any(abs(colSums(abs(x)) - 1) > 10e-3) && (perc == T)){
        warning("turning the sum of loadings to percentage is not appropriate")
        perc = FALSE
        digits = 3
      }
      out = apply(x, 2, function(y, ii) tapply(y, ii, sum), ii = ind)
    if (only.nonzero)
      out = out[rowSums(abs(out)) > 10e-5, ]
    }
    }
    if(is.numeric(digits)){
      if(isTRUE(perc))
        out = round(100*out, digits)
      else
        out = round(out, digits)
      }
    else{      
      if(isTRUE(perc))
        out = 100*out
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
#' @param object An spca object.
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
#' @param ... Further arguments; currently ignored.
#' 
#' @return If rtn = TRUE, a numerical matrix with the summaries.
#' @seealso Examples in \code{\link{lsspca}}.
#' @export
#' @method summary spca
summary.spca <- function(object, cols, contribution = TRUE, 
                         variance_metrics = "cumulative_relative",
                         minload = FALSE, rtn = FALSE, prn = TRUE, 
                         thrsehcard = 0.001, ...) {
  x = object
  # Validation
  if (!is.spca(x)) {
    stop("summary.spca works only for spca objects")
  }
  
  if (is.vector(x$loadings)) {
    x$loadings <- as.matrix(x$loadings)
  }
  
  # Determine columns
  if (missing(cols)) {
    cols <- 1:min(ncol(x$loadings), length(x$vexpPC))
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
    VEXP = x$vexp * 100,
    CVEXP = cumsum(x$vexp) * 100
  )
  
  # Add variance comparison metrics based on user choice
  if (variance_metrics %in% c("relative", "both")) {
    out <- rbind(out, RVEXP = x$vexp / x$vexpPC * 100)
  }
  
  if (variance_metrics %in% c("cumulative_relative", "both")) {
    out <- rbind(out, RCVEXP = cumsum(x$vexp) / cumsum(x$vexpPC) * 100)
  }
  
  # Add cardinality
  out <- rbind(out, Card = apply(abs(x$loadings) > thrsehcard, 2, sum))
  
  # Add convergence if available
  if (!is.null(x$conv)) {
    out <- rbind(out, Converged = x$conv)
  }
  
  # Add minimum loading/contribution if requested
  if (minload) {
    if (contribution) {
      min_vals <- get.minload(x, perc = TRUE) * 100
      min_label <- "MinCont"
    } else {
      min_vals <- get.minload(x, perc = FALSE)
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
