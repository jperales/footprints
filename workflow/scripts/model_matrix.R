# point of this file:
# - use the zscores to create a linear model
library(dplyr)
library(broom)
library(modules)
b = import('../../util/ebits/base')
io = import('../../util/ebits/io')
ar = import('../../util/ebits/array')


#' Fits a linear model on Z-scores
#'
#' @param zdata  A list with the zscore matrix and index object
#' @return       The coefficients matrix [gene x pathway]
zscore2model = function(zdata, hpc_args=NULL) {
    index = zdata$index
    zscores = t(zdata$zscores) * index$sign
    temp=apply(is.na(zscores),2,sum)
    temp=temp[temp<10]
    zscores=zscores[,names(temp)]
    # fit model to pathway perturbations
    pathwayX = t(ar$mask(index$pathway)) + 0
    # NOTE: add EGFR>PI3K link here?
    #pathwayX["PI3K",] = pathwayX["EGFR",] + pathwayX["PI3K",]
    pathwayX["MAPK",] = pathwayX["EGFR",] + pathwayX["MAPK",]
    pathwayX["NFkB",] = pathwayX["TNFa",] + pathwayX["NFkB",]
    pathwayX=t(pathwayX)
    mod = lm(zscores ~ 0 + pathwayX) %>% broom::tidy() %>%
      transmute(gene = response,
                pathway = sub("^pathwayX", "", term),
                zscore = estimate,
                p.value = p.value) %>%
      mutate(adj.p = p.adjust(p.value, method="fdr"))


    #NOTE: Original shuttle lm replaced by stats::lm and broom::tidy() reformat
    #---
    #index = zdata$index
    #zscores = t(zdata$zscores) * index$sign

    # fit model to pathway perturbations
    #pathway = t(ar$mask(index$pathway)) + 0
    #pathway["EGFR",] = pathway["EGFR",] + pathway["MAPK",] + pathway["PI3K",]
    #pathway["TNFa",] = pathway["TNFa",] + pathway["NFkB",]

    #mod = st$lm(zscores ~ 0 + pathway, data=index, min_pts=30, atomic="pathway",
    #            hpc_args=hpc_args) %>%
    #    transmute(gene = zscores,
    #              pathway = sub("^pathway", "", term),
    #              zscore = estimate,
    #              p.value = p.value) %>%
    #    mutate(adj.p = p.adjust(p.value, method="fdr"))
    #---

    zfit = ar$construct(zscore ~ gene + pathway, data=mod)
    pval = ar$construct(p.value ~ gene + pathway, data=mod)

    # filter zfit to only include top 100 genes per pathway
    model = zfit
    model[apply(pval, 2, function(p) !b$min_mask(p, 100))] = 0

    list(assocs=mod, model=model, zfit=zfit)
}

if (is.null(module_name())) {
    ZDATA = commandArgs(TRUE)[1] %or% "../../data/zscores.RData"
    OUTFILE = commandArgs(TRUE)[2] %or% "../../model/model_matrix"

    # load speed data, index; filter for train set only
    zdata = io$load(ZDATA)
    result = zscore2model(zdata, hpc_args=NULL)
    # save resulting object
    save(result, file=paste0(OUTFILE,'.RData'))
    write.csv(result$model,paste0(OUTFILE,'.csv'))
    write.csv(result$zfit,paste0(OUTFILE,'_full','.csv'))
}
