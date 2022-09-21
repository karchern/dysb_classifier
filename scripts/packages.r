library(hms)
library(tidyverse)
library(readxl)
library(ggplot2)
library(SIAMCAT)
library(data.table)
library(vegan)
library(RColorBrewer)

parseArgs <- function() {

    ## Collect arguments
    args <- commandArgs(TRUE)

    ## Parse arguments (we expect the form --arg=value)
    parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
    argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
    argsL <- as.list(as.character(argsDF$V2))
    names(argsL) <- argsDF$V1

    return(argsL)

}

scale_color_shape_combined <- function(ggplotObject, var, varName, numShapes = round(length(levels(var))/3, 0), colors = c("#7fc97f", "#beaed4", "#fdc086")) {
    if (!is.factor(var)) {
        print("Supplied variable needs to be a factor!")
        exit()
    }
    
#     #print(numColors)
#     print(numShapes)
#     print(varName)
#     print(head(var))
#     print(head(levels(var)))
#     print(length(levels(var)))
    
#     print("####")
#     print("####")
#     print(rep(colors, length.out = length(levels(var))))
#     print(rep(1:numShapes, length.out = length(levels(var))))
    
    ggplotObject + 
           scale_colour_manual(name = varName,
                       labels = levels(var),
                       values = rep(colors, length.out = length(levels(var)))) +
    scale_shape_manual(name = varName,
                   labels = levels(var),
                   values = rep((1:numShapes), length.out = length(levels(var)))) %>%
    return()
}

getFeatures <- function(dataMetaAdaptedTruncWideAll, metaDataWGS) {

    # Overall Richness
    richness <- apply(dataMetaAdaptedTruncWideAll, 2, function(x) sum(x>0))
    richnessWithCutoff <- apply(dataMetaAdaptedTruncWideAll, 2, function(x) sum(x>1E-5))
                      
    # what is less prevalent than 5% in healthy controls? Sum up prevalences and add as feature
    HC <- dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% (metaDataWGS %>% 
                                                                                          filter(caseControls == "control") %>% 
                                                                                          mutate(r = rownames(.)) %>% 
                                                                                          pull(r))]
    absentInHC <- apply(HC, 1, function(x) mean(x>0) < 0.05)
    absentInHCWithCutoff <- apply(HC, 1, function(x) mean(x>1E-5) < 0.05)                    
    richnessAbsentInHC <- apply(dataMetaAdaptedTruncWideAll[absentInHC, ], 2, function(x) sum(x>0))
    richnessWithCutoffAbsentInHCWithCutoff <- apply(dataMetaAdaptedTruncWideAll[absentInHCWithCutoff, ], 2, function(x) sum(x>1E-5))                            
    # what is less prevalent than 5% in healthy controls? Sum up abundances and add as feature
    cumAbundanceAbsentInHC <- apply(dataMetaAdaptedTruncWideAll[absentInHC, ], 2, function(x) sum(x))                              
    cumAbundanceAbsentInHCWithCutoff <- apply(dataMetaAdaptedTruncWideAll[absentInHCWithCutoff, ], 2, function(x) sum(x))                                                                  
    # what is more prevalent than 50% in healthy controls? Sum up their prevalences and add as feature
    presentInHC <- apply(HC, 1, function(x) mean(x>0) > 0.5) 
    presentInHCWithCutoff <- apply(HC, 1, function(x) mean(x>1E-5) > 0.5)                          
    richnessPresentInHC <- apply(dataMetaAdaptedTruncWideAll[presentInHC, ], 2, function(x) sum(x>0))
    richnessWithCutoffPresentInHCWithCutoff <- apply(dataMetaAdaptedTruncWideAll[presentInHCWithCutoff, ], 2, function(x) sum(x>1E-5))                                 
    # what is more prevalent than 80% in healthy controls? Sum up their abundances and add as feature
    cumAbundancePresentInHC <- apply(dataMetaAdaptedTruncWideAll[presentInHC, ], 2, function(x) sum(x))
    cumAbundancePresentInHCWithCutoff <- apply(dataMetaAdaptedTruncWideAll[presentInHCWithCutoff, ], 2, function(x) sum(x))
                                               
    # relAbundance of top taxon
    topTaxonRelAb <- apply(dataMetaAdaptedTruncWideAll, 2, function(x) max(x))   
    # meanRelAbundance of 5% of non-zero taxa (skew of abundance distribution)
    topTaxaRelAb <- apply(dataMetaAdaptedTruncWideAll, 2, function(x) {
        x <- x[x!=0]
        l <- length(x)
        return(mean(x[as.integer(l*0.95):l]))            
    }) 
                                                   

    morePredictorNames <- c('richness',
                            "richnessWithCutoff",
                           'richnessAbsentInHC',
                           'richnessWithCutoffAbsentInHCWithCutoff',
                           'cumAbundanceAbsentInHC',
                           'cumAbundanceAbsentInHCWithCutoff',
                           'richnessPresentInHC',
                           'richnessWithCutoffPresentInHCWithCutoff',
                           'cumAbundancePresentInHC',
                           'cumAbundancePresentInHCWithCutoff',
                           'topTaxonRelAb',
                           'topTaxaRelAb')                                 

    allDat <-   list(richness,
                  richnessWithCutoff,
                  richnessAbsentInHC,
                  richnessWithCutoffAbsentInHCWithCutoff,   
                  cumAbundanceAbsentInHC,
                  cumAbundanceAbsentInHCWithCutoff,
                  richnessPresentInHC,
                  richnessWithCutoffPresentInHCWithCutoff,   
                  cumAbundancePresentInHC,
                  cumAbundancePresentInHCWithCutoff,
                  topTaxonRelAb,
                  topTaxaRelAb)   
    names(allDat) <- morePredictorNames                                 

    return <- list(allDat, morePredictorNames)

}