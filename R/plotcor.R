## 
mkf = function(x){
  g = x[1]
  levs = rep(0, length(x))
  labs = g
  k = 1
  for(i in 1:length(x)){
    if(x[i] == g)
      levs[i] = k
    else{
      k = k + 1
      g = x[i]
      levs[i] = k
      labs = c(labs, x[i])
    }
  }
  factor(levs, labels = labs)
}

#' Correlation heatmap with optional grouping and annotations
#'
#' Draws a correlation heatmap from a square matrix `S`, optionally reordering variables
#' by groups, adding block separators (and subgroup rectangles), and overlaying numeric correlation labels.
#'
#' @param S A square numeric matrix (typically a correlation matrix) with values in \eqn{[-1, 1]}. Rows/columns are reordered if `groups` is provided.
#' @param axis_labels Controls axis labels. `TRUE` uses column/row names if available (otherwise uses `V1, \dots, Vp`); `FALSE` or `"none"` removes labels. A character/factor vector can be supplied to use custom labels.
#' @param groups Optional grouping information used to reorder variables and (optionally) draw separators. Accepted inputs: `FALSE` (no grouping), a permutation vector of indices `1:p`, a character/factor vector of length `p` (group labels), or a list of
#'   integer vectors (one per group).
#' @param subgroups Optional subgroup specification (same accepted formats as `groups`). When provided and groups are separated, subgroup rectangles are drawn.
#' @param separate_groups Logical or numeric. If positive, draws separators between groups with the specified line width.
#' @param group_sepline Color for group separator lines.
#' @param bnw Logical. If `TRUE`, uses a black/white fill scale; otherwise uses a diverging scale.
#' @param addlabels Logical. If `TRUE`, overlays numeric correlation labels inside cells.
#' @param thresholdFontColLabels Numeric. Threshold used to choose black/white text for overlaid labels (relative to the displayed values).
#' @param label_size Numeric. Size of overlaid labels.
#' @param label_face Character. Font face for overlaid labels (e.g., `"plain"`, `"bold"`).
#' @param lowcol,highcol Length-2 character vectors. Colors used for the diverging fill scale (only the first element is used).
#' @param axis_lab_size Numeric. Axis label text size.
#' @param bottom_axis_labels,side_axis_labels Logical. Whether to show x- and y-axis labels.
#' @param add_group_names Logical. If `TRUE`, adds group names above the plot (requires `separate_groups > 0` and grouped input).
#' @param group_names Optional character vector of group names (length must match the number of groups). If `NULL` and groups are named, names are used.
#' @param group_names_size Numeric. Size of group name labels.
#' @param group_names_face Character. Font face for group name labels.
#' @param group_names_vjust Numeric. Vertical adjustment for group name labels.
#' @param expandTop Numeric vector. Expansion applied to the top of the plot when adding group names (passed to `scale_y_discrete(expand = expansion(...))`).
#' @param clipPlotOff Logical. If `TRUE`, sets `coord_cartesian(clip = "off")` to allow group names to extend beyond the panel.
#' @param subgroup_sepline_colour Color for subgroup rectangle outlines.
#' @param subgroup_sepline_size Numeric. Line width for subgroup rectangle outlines.
#' @param plt Logical. If `TRUE`, prints the plot.
#' @param rtn_plot Logical. If `TRUE`, returns the `ggplot` object.
#'
#' @return If `rtn_plot = TRUE`, a `ggplot` object; otherwise `NULL` (invisibly). If`plt = TRUE`, the plot is printed.
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[reshape2]{melt}}
#'
#' @export
plotcor <- function (S, axis_labels = TRUE, groups = FALSE, 
                      subgroups = FALSE, 
                     separate_groups = (!(groups[[1]] == FALSE))[1], 
                     group_sepline = "black", bnw = FALSE, 
                     addlabels = FALSE, thresholdFontColLabels = 0.33, 
                     label_size = 5, label_face = "plain",  
                     lowcol = c("firebrick", "darkgreen"), 
                     highcol = c("dodgerblue3", "darkblue"),
                     axis_lab_size = 10, 
                     bottom_axis_labels = FALSE,
                     side_axis_labels = TRUE,
                     add_group_names = F, group_names = NULL,
                     group_names_size = 4, group_names_face = "plain",
                     group_names_vjust = 0, 
                     expandTop = c(0, 0.1),
                     clipPlotOff = FALSE,
                     subgroup_sepline_colour  = "black", 
                     subgroup_sepline_size = 1, 
                     plt = TRUE, 
                     rtn_plot = FALSE
                     ) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
    stop("ggplot2 package needed, run install(ggplot2)")
#  else require("ggplot2")
  if (!requireNamespace("reshape2", quietly = TRUE))
    stop("reshape2 package needed, run install(reshape2)")
# #  else require("reshape2")
  
  p = ncol(S)
  
  #if (is.null(rownames(S) && is.null(colnames(S)
                                     
  if((!(groups[[1]][1] == FALSE)) &  (!is.list(groups))){
    if(is.vector(groups))
      groups = mkf(groups)
    if(is.factor(groups)){
      nam = levels(groups)
      groups = lapply(levels(groups), 
                             function(m, x) which(x == m), x = groups)
    names(groups) = nam
      }
    } 
  if( (!(subgroups[[1]][1] == FALSE)) &  (!is.list(subgroups))){
  if(is.vector(subgroups))
    subgroups = mkf(subgroups)
  if(is.factor(subgroups))
    subgroups = lapply(levels(subgroups), 
                           function(m, x) which(x == m), x = subgroups)
}  
  
  
  if(is.list(groups))
    ind_list = groups
  else 
    ind_list = NA
  
  if (groups[[1]][1] == FALSE){
    ind = 1:p
    separate_groups = FALSE
  }
  else 
      if (is.vector(groups) & !is.list(groups)){
        ind = groups
        separate_groups = FALSE
        if (!all(sort(ind) == 1:p))
          stop("There are missing indices groups")
      }
      else{
        if (!is.list(groups) && length(groups) < p)
      groups = list(groups, (1:p)[-groups])
    if (is.list(groups)){
      nblocks = length(groups)
#separates blocks==================================      
      if (separate_groups == TRUE){
        if(add_group_names){
          if((is.list(groups)) && !is.null(names(groups)) &&
             is.null(group_names))
            group_names = names(groups)
          else
            if(is.factor(groups) &&
               is.null(group_names))  group_names = levels(groups)
            else
              if (is.vector(groups)&&
                  is.null(group_names))  group_names = unique(groups)
        }
        nvar = sapply(groups, length) #c(7, 6, 7)  
        cv = c(0,cumsum(nvar))
        nv = length(nvar) + 1
        xM = cv[-nv] + 0.5
        yM = p - cv[-nv] + 0.5
        xm = cv[-1] + 0.5
        ym = p - cv[-1] + 0.5 
        xmm = sapply(ind_list[-nblocks], max) + 0.5
      }##
      ind = unlist(groups)
      if (!all(sort(ind) == 1:p))
        stop("There are missing indices groups")
    }
  }
  if (!is.vector(ind))
    stop("something wrong with vargrouop_list")
    else
      S = S[ind, ind]

  if(axis_labels == FALSE){
    bottom_axis_labels = FALSE  
    side_axis_labels = FALSE  
  }
  
##  axislables = TRUE
  if ( (axis_labels[1] == FALSE) || (axis_labels[1] == "none")) {
    #dimnames(S) = NULL qui problema con melt sotto
  }
  else{
    if (is.vector(axis_labels) & (length(axis_labels) > 1)) {
      if (length(axis_labels) < p){
        print("fewer labels than variables, axis_labels will be V1,..,Vp")
      }
    }
    else 
      if (is.factor(axis_labels)) 
        axis_labels <- as.character(axis_labels)
    if (axis_labels[1] == TRUE){ 
        if(all(is.null(dimnames(S))))
          axislables = paste0("V", 1:p)
        else{
          if(is.null(colnames(S)))
            axis_labels = paste0("V", 1:p)
            axis_labels = paste0("V", 1:p)
        }
  }
} 
  
## plot 
  dat <- reshape2::melt(S[, p:1])
  q <- ggplot2::ggplot(data = dat, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_raster() + ggplot2::xlab( "") + ggplot2::ylab( "") +
    ggplot2::theme(axis.text = ggplot2::element_text(size = axis_lab_size)) #element_text(angle=90, hjust=1)
  if (bottom_axis_labels == FALSE)
    q = q +  ggplot2::theme(axis.text.x = element_blank(), axis.ticks = element_blank()) 
  if (side_axis_labels == FALSE)
    q = q +  ggplot2::theme(axis.text.y = element_blank())
  
  if (bnw == TRUE){
    q = q + ggplot2::theme_bw()  + 
      ggplot2::scale_fill_gradient(limits = c(-1,1), high = "black", low = "white", guide =  ggplot2::guide_legend(title ="Cor")) +
      ggplot2::theme(legend.position="bottom",
            legend.key = ggplot2::element_rect(color = "black", size = 1, linetype = 1))
  }
  else 
    if (!(is.na(lowcol[1]) | is.na(highcol[1])))    
      q = q + ggplot2::scale_fill_gradient2(limits = c(-1, 1), low = lowcol[1], high = highcol[1], guide = ggplot2::guide_legend(title ="Cor"))  + 
    ggplot2::theme(legend.position="bottom",
                   legend.key = ggplot2::element_rect(color = "black", size = 1, linetype = 1))
  else
    q = q + ggplot2::scale_fill_gradient2(limits = c(-1, 1), 
                      guide =  ggplot2::guide_legend(title ="Cor"))  + 
    ggplot2::theme(legend.position="bottom",
                   legend.key = ggplot2::element_rect(color = "black", size = 1, linetype = 1))
  
  if (separate_groups > 0){
    hline_df = data.frame(y = yM[-1], yend = yM[-1], x = 0.5, xend = p +0.5) 
    q = q + ggplot2::geom_segment(inherit.aes = FALSE, data = hline_df, aes(x = x, xend = xend, y = y, yend = yend),  linewidth = separate_groups)
    
    nvlines = length(xmm)
    vline_df = data.frame(x = xmm, y = rep(0.5, nvlines), yend = rep(p, nvlines) + 0.5)
    
    q <- q +  ggplot2::geom_segment(inherit.aes = FALSE, data = vline_df, aes(x = x, xend = x, y = y, yend = yend),  colour= group_sepline, lwd = separate_groups)
    
    if (subgroups[[1]][1] != FALSE){
      if (!is.list(subgroups) && length(subgroups) < p)
        subgroups = list(subgroups, (1:p)[-subgroups])
      if (is.list(subgroups)){
        xm = sapply(subgroups, min) - 0.5
        xM = sapply(subgroups, max) + 0.5
        ym = p - sapply(subgroups, max) + 0.5
        yM = p - (sapply(subgroups, min) ) + 1.5
        
        
        
        q = q +  ggplot2::annotate(geom = "rect", xmin = xm, xmax = xM, ymin = ym, ymax = yM, colour = subgroup_sepline_colour, size = subgroup_sepline_size, fill="transparent")
      } 
    }
    if(separate_groups > 0){
      if(add_group_names){
        if(length(group_names) != nblocks){
          message("cannot add block names, too few or too many")
          group_names = NA
        }
        
        else{ 
          if(add_group_names){
            st = c(0, cumsum(nvar)[-nblocks]) + 0.5
            hjustb = 0
            if (is.list(ind_list)){
              st = sapply(ind_list, mean)
              hjustb = 0.5
            }
            top_y <- utils::tail(levels(dat$Var2), 1)
            
            q <- q +
              ggplot2::annotate(
                "text",
                x = st, y = top_y,
                label = group_names,
                hjust = hjustb,
                vjust = -0.6 + group_names_vjust,
                size = group_names_size,
                fontface = group_names_face
              )
            if (expandTop[1] > 0)
              q = q +  scale_y_discrete(expand = expansion(mult = c(0.01, expandTop)))
            
            
            if(clipPlotOff) q <- q + ggplot2::coord_cartesian(clip = "off")
          }
        }
      }
    }
  }
  if (bottom_axis_labels == FALSE)
    q = q +  ggplot2::theme(axis.text.x=element_blank())
  else
    if (axis_labels[1] != FALSE) 
      q = q + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  if (side_axis_labels == FALSE)
    q = q +  ggplot2::theme(axis.text.y = element_blank())
  
  if (addlabels == TRUE){
    myf = function(x) ifelse (rev(x) < thresholdFontColLabels, "black", "white")
    mic = apply(t(S), 2, myf)
    myl = function(x) round(rev(x), 1)
    labs =  apply(S, 1, myl)
    q = q + ggplot2::geom_text(aes(x = rep(1:p, each = p), y = rep(1:p, p), label = c(labs)),  colour = c(mic), size = label_size, fontface = label_face)
  }
  q = q + ggplot2::theme(panel.background = element_rect(fill = "white", 
                                                colour = "black", linetype = 1))
  if (plt == TRUE)
    print(q)

  if (rtn_plot == TRUE)
    return(q)
}
