library(hms)
library(tidyverse)
library(readxl)
library(ggplot2)
library(SIAMCAT)

load("/g/scb2/zeller/karcher/dysb_classif/data/siamcat/profiles_merged_with_metadata.rimage")

# Low-prevalence filtering
filterBool <- apply(dataMetaAdaptedTruncWideAll, 1, function(x) mean(x>0) > 0.05)
print(table(filterBool))                    
dataMetaAdaptedTruncWideAll <- dataMetaAdaptedTruncWideAll[filterBool, ]
stopifnot(all(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS)))                    

args<-commandArgs(TRUE)
cond <- args[1]

    # Train on everything but a given study
    meta.train <- metaDataWGS %>% 
        # For training include cohort studies/only-control studies
        #filter(condition != cond | condition == "No condition") %>%
	# Actually play around by removing cohort studies (for now)
	filter(condition != cond & condition != "No condition" & condition != "T2D") %>%
        as.data.frame()
    
    block <- 'condition'

    label <- create.label(meta=meta.train ,
            label='caseControls', case='case')
    
    if (! all(rownames(meta.train) %in% colnames(dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)])) || ! all(colnames(dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)]) %in% rownames(meta.train))){
        print("No perfect overlap between training metadata and profiles. Exiting")
        sadsadsa
    }
    
    # create SIAMCAT object
    sc.obj.train <- siamcat(feat=dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)], meta=meta.train,
                            label=label)
    
    # normalize features
    sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
        norm.param=list(log.n0=1e-05, sd.min.q=0),feature.type = 'original')
    # Create data split
    sc.obj.train <- create.data.split(sc.obj.train,
        num.folds = 5, num.resample = 1, inseparable = block)
    # train LASSO model
    sc.obj.train <- train.model(sc.obj.train, method='randomForest')

    # apply trained models to held-out dataset
    meta.test <- metaDataWGS %>% 
        filter(condition==cond) %>%
        as.data.frame()
    sc.obj.test <- siamcat(feat=dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.test)], meta=meta.test,
                            label='caseControls', case='case')

    # make holdout predictions
    sc.obj.test <- make.predictions(sc.obj.train, 
                                    siamcat.holdout = sc.obj.test,
                                    normalize.holdout = TRUE)
    sc.obj.test <- evaluate.predictions(sc.obj.test)
    
    sc.obj.train.LOCO <- sc.obj.train

    #######################
    # within-condition CV #
    #######################

    meta.train <- metaDataWGS %>%
        # For training include cohort studies/only-control studies
        #filter(condition != cond | condition == "No condition") %>%
        # Actually play around by removing cohort studies (for now)
        filter(condition == cond) %>%
        as.data.frame()
    
    #block <- 'condition'

    label <- create.label(meta=meta.train ,
            label='caseControls', case='case')
    
    if (! all(rownames(meta.train) %in% colnames(dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)])) || ! all(colnames(dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)]) %in% rownames(meta.train))){
        print("No perfect overlap between training metadata and profiles. Exiting")
        sadsadsa
    }
    
    # create SIAMCAT object
    sc.obj.train <- siamcat(feat=dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)], meta=meta.train,
                            label=label)
    
    # normalize features
    sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
        norm.param=list(log.n0=1e-05, sd.min.q=0),feature.type = 'original')
    # Create data split
    sc.obj.train <- create.data.split(sc.obj.train,
        num.folds = 5, num.resample = 1)
    # train LASSO model
    sc.obj.train <- train.model(sc.obj.train, method='randomForest')
    sc.obj.train <- make.predictions(sc.obj.train)
    sc.obj.train <- evaluate.predictions(sc.obj.train)

save(sc.obj.train.LOCO, sc.obj.test, sc.obj.train, cond, file =  str_c("/g/scb2/zeller/karcher/rtmp/", cond))
