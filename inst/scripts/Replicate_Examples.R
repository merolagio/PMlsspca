
#MSSCQ data ===========================
load( "Data/msscq.Rdata", verbose = T)
# Loading objects:
# ms
# ms_lookup
# ms_scalesh_fac

ms_scalesh_names = unique(ms_lookup$scaleSH)

ms_n = nrow(ms)
ms_p = ncol(ms)

ms_r = cor(ms)

plotcor(ms_r, groups = ms_scalesh_fac,   separate_groups = T, axis_labels = F, add_group_names = T, expandTop = 0.05, group_names_vjust = -0.05)
##scree and qq plots===========

msr_ee = eigen(ms_r, symmetric = T)
screeplot(msr_ee$values)
wachterqq(msr_ee$values, p = ms_p, n = ms_n, cor = T, nfit_line = 96)


##PCA==============

ms_pca = pca(ms, ncomps = 4)
summary(ms_pca)

# plot contributions
ms_pca_contr_pl = plot.spca(ms_pca, returnplot = T, vargroups = ms_scalesh_fac, varnames = F, colourscale = "ggplot", legendPosition = "right", stripnames = paste("Component", 1:4), produceplot = F)

ms_pca_contr_pl

##SPCA===================

#slow 17k obs x 100 vars hard for leaps::regsubsets

ms_pspcas95 = lsspca(ms, alpha = 0.95, maxcard = 0, ncomps = 4,  method = "p", varselection = "s", scalex = FALSE) 

summary(ms_pspcas95, variance_metrics = "both")

ms_pspcas95_load_pl = plot.spca(ms_pspcas95, nplot = 4, onlynonzero = T, varnames = F, legendPosition = "right", vargroups = ms_scalesh_fac, colourscale = "ggplot", returnplot = T, produceplot = F)  +  labs(x = "item")

ms_pspcas95_load_pl 


# correlation table
tab_cors = round(rbind(ms_pspcas95$corComp, diag(cor(ms_pca$scores, ms_pspcas95$scores))), 2)
tab_cors

#CRIME=================

#load("Data/cr.RData", verbose = T)

data("cr", package = "PMlsspca")

cr_p = ncol(cr)
cr_n = nrow(cr)

cr_r = cor(cr)

#scree and qq plots===========

##CLUSTER======================
cr_hc <- hclust(as.dist(1 - abs(cr_r)), method = "average")
cr_hcord <- cr_hc$order

crhc = cr[, cr_hcord]
crhc_r = cor(crhc)
colnames(crhc)[1:10]
crhc_cor_pl = plotcor(crhc_r, rtn = T, axis_labels = F)


cr_ee = eigen(cr_r)

screeplot(cr_ee$values)
wachterqq(cr_ee$values, p = ms_p, n = ms_n, cor = T, nfit_line = -3)

##PCA===============

cr_pca = pca(cr, ncomps = 4)
summary(cr_pca)

##SPCA===============

cr_pspcas95 = lsspca(cr, alpha = 0.95, ncomps = 4,  method = "p", varselection = "s", scalex = FALSE) 

##tables====================

###           summaries 
su = summary(cr_pspcas95, variance_metrics = "both", rtn = T, prn = F)
PCs = c(diag(cor(cr_pca$scores, cr_pspcas95$scores)))
names(x) = colnames(su)
round(rbind(su, PCs), 2)

###           contributions 

tab_contr = cr_pspcas95$contributions[, 1:2]
ind = rowSums(tab_contr) != 0
tab_contr = tab_contr[ind, ]
df_contr = data.frame(cr_lookup[ind, c(1, 3)], as.data.frame(round(100*tab_contr, 0)))
ind = c(which(tab_contr[, 1] != 0), which(tab_contr[, 1] == 0))
df_contr =  df_contr[ind, ]
df_contr
