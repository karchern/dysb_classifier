source('/g/scb2/zeller/karcher/dysb_classif/scripts/packages.r')
args <- parseArgs()
load(str_c(args[['BASEDIR']], "data", "profiles_merged_with_metadata.rimage", sep = "/"))

# Low-prevalence filtering
if (args[['lowAbFiltering']]) {
    filterBool <- apply(dataMetaAdaptedTruncWideAll, 1, function(x) mean(x>0) > 0.05)
    dataMetaAdaptedTruncWideAll <- dataMetaAdaptedTruncWideAll[filterBool, ]
    stopifnot(all(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS)))   
}
                        
if (args[['subsample']]) {
    set.seed(1337)
    metaDataWGS <- metaDataWGS %>% 
        mutate(sampleID = rownames(.)) %>%
        group_by(condition) %>% nest() %>% 
        mutate(data = map(data, function(x) {

    if (dim(x)[1] > 200) {
        x %>%
        sample_n(200, replace = F) %>%
        return()
    } else {
        x %>% return()
    }
    })) %>%
    unnest() %>%
    as.data.frame()
    rownames(metaDataWGS) <- metaDataWGS$sampleID
    metaDataWGS$sampleID <- NULL

    dataMetaAdaptedTruncWideAll <- dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(metaDataWGS)]
    dataMetaAdaptedTruncWideAll  <- dataMetaAdaptedTruncWideAll[, match(rownames(metaDataWGS), colnames(dataMetaAdaptedTruncWideAll))]
    #print(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS))
    #print(dataMetaAdaptedTruncWideAll %>% dim())
    #print(metaDataWGS %>% dim())
    stopifnot(all(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS)))
}
                        
if (args[['genusLevel']]) {
    tax <- read_tsv('../data/motus2.6_taxonomy_NCBI.tax') %>%
    select(mOTUs_ID, Genus, profiled)
    tax <- tax %>% 
        mutate(profile_corrected = map_chr(profiled, function(x) {
        str_c(str_split(x, " ")[[1]][2:length(str_split(x, " ")[[1]])], collapse = " ") %>%
            return()
        })) 

    # TODO Figure out!
    # We loose a few dozen mOTUs (< 0.02% of all mOTUs) cause they're not there in the tax file, not sure why...
    tmp <- dataMetaAdaptedTruncWideAll %>% 
    mutate(profile_corrected = rownames(.)) %>% 
    inner_join(tax %>% select(profile_corrected, Genus), by = 'profile_corrected') %>% 
    select(-profile_corrected) %>% 
    pivot_longer(-Genus) %>% 
    rename(sampleID = name, relAb = value)

    tmp2 <- tmp %>% group_by(sampleID, Genus) %>%
    summarize(relAb = sum(relAb)) %>%
    pivot_wider(id_cols = Genus, names_from = sampleID, values_from = relAb)

    tmp3 <- tmp2 %>%
    as.data.frame()
    rownames(tmp3) <- tmp3$Genus
    tmp3$Genus <- NULL
    dataMetaAdaptedTruncWideAll <- tmp3   
    dataMetaAdaptedTruncWideAll <- dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(metaDataWGS)]
    dataMetaAdaptedTruncWideAll  <- dataMetaAdaptedTruncWideAll[, match(rownames(metaDataWGS), colnames(dataMetaAdaptedTruncWideAll))]
    #print(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS))
    #print(dataMetaAdaptedTruncWideAll %>% dim())
    #print(metaDataWGS %>% dim())
    stopifnot(all(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS)))
}
                        
if (args[['featureEng']]) {
allDat <- getFeatures(dataMetaAdaptedTruncWideAll, metaDataWGS) 
morePredictorNames <- allDat[[2]]
allDat <- allDat[[1]]    
nn <- names(allDat)                           
stopifnot(all(map_dbl(allDat, length) == dim(dataMetaAdaptedTruncWideAll)[2]))    
    
tmpData <- matrix(unlist(allDat), byrow = T, nrow = length(allDat)) %>%
                                 as.data.frame()
colnames(tmpData) <- colnames(dataMetaAdaptedTruncWideAll)
rownames(tmpData) <- nn   
tmpData <- t(tmpData)
stopifnot(all(rownames(tmpData) == rownames(metaDataWGS)))
metaDataWGS <- metaDataWGS %>% 
                                 as.data.frame() %>% 
                                 mutate(sampleID = rownames(.)) %>% 
                                 left_join(tmpData %>% 
                                           as.data.frame() %>% 
                                           mutate(sampleID = rownames(.)), 
                                           by ='sampleID')                                
rownames(metaDataWGS) <- metaDataWGS$sampleID
metaDataWGS$sampleID <- NULL                                 
}                       
                 
#args<-commandArgs(TRUE)
cond <- args[['condition']]

# Train on everything but a given study
meta.train <- metaDataWGS %>% 
    # For training include cohort studies/only-control studies
    #filter(condition != cond | condition == "No condition") %>%
    # Actually play around by removing cohort studies (for now)
    filter(condition != cond & condition != "No condition" & condition != "T2D") %>%
    as.data.frame()
                 


label <- create.label(meta=meta.train ,
        label='caseControls', case='case')
save.image('/g/scb2/zeller/karcher/tmp.rsave')
asdssa

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
                                 
if (args[['block']] == 'condition') {
    sc.obj.train <- create.data.split(sc.obj.train,
    num.folds = 5, num.resample = 1, inseparable = 'condition')
} else {
    sc.obj.train <- create.data.split(sc.obj.train,
    num.folds = 5, num.resample = 1)
}                        

# potentially: Add features
if (args[['featureEng']]) {
    sc.obj.train <- add.meta.pred(sc.obj.train, morePredictorNames)
}                                 
                        
# train model
if (args[['MLalgorithm']] == "RF") {
  sc.obj.train <- train.model(sc.obj.train, method='randomForest')  
} else if (args[['MLalgorithm']] == "lasso") {
  sc.obj.train <- train.model(sc.obj.train, method='lasso')  
} else if (args[['MLalgorithm']] == "lasso_ll") {
  sc.obj.train <- train.model(sc.obj.train, method='lasso_ll')  
} else {
  print("ML Algorithm undefined. exciting.")
  stop()
}

# apply trained models to held-out dataset
meta.test <- metaDataWGS %>% 
    filter(condition==cond) %>%
    as.data.frame()
                        
sc.obj.test <- siamcat(feat=dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.test)], meta=meta.test,
                        label='caseControls', case='case')


# Explicitly normalize features in test set so that you can then add metadata features to test set                                
sc.obj.test <- normalize.features(sc.obj.test, norm.method = 'log.std',
    norm.param=norm_params(sc.obj.train),feature.type = 'original')      
                                 
# Add features
if (args[['featureEng']]) {
    sc.obj.test <- add.meta.pred(sc.obj.test, morePredictorNames)
}   
                                                          
# make holdout predictions. DONT NORMALIZE AGAIN
sc.obj.test <- make.predictions(sc.obj.train, 
                                siamcat.holdout = sc.obj.test,
                                normalize.holdout = FALSE)                                 
                                 
# if (args[['featureEng']]) {
#     sc.obj.test <- add.meta.pred(sc.obj.test, morePredictorNames)
# }   
# # make holdout predictions
# sc.obj.test <- make.predictions(sc.obj.train, 
#                                 siamcat.holdout = sc.obj.test,
#                                 normalize.holdout = TRUE)

                                                                  
sc.obj.test <- evaluate.predictions(sc.obj.test)

# #######################
# # within-condition CV #
# #######################

# meta.train <- metaDataWGS %>%
#     # For training include cohort studies/only-control studies
#     #filter(condition != cond | condition == "No condition") %>%
#     # Actually play around by removing cohort studies (for now)
#     filter(condition == cond) %>%
#     as.data.frame()

# #block <- 'condition'

# label <- create.label(meta=meta.train ,
#         label='caseControls', case='case')

# if (! all(rownames(meta.train) %in% colnames(dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)])) || ! all(colnames(dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)]) %in% rownames(meta.train))){
#     print("No perfect overlap between training metadata and profiles. Exiting")
#     sadsadsa
# }

# # create SIAMCAT object
# sc.obj.train <- siamcat(feat=dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.train)], meta=meta.train,
#                         label=label)

# # normalize features
# sc.obj.train <- normalize.features(sc.obj.train, norm.method = 'log.std',
#     norm.param=list(log.n0=1e-05, sd.min.q=0),feature.type = 'original')
# # Create data split
# sc.obj.train <- create.data.split(sc.obj.train,
#     num.folds = 5, num.resample = 1)
# # train LASSO model
# sc.obj.train <- train.model(sc.obj.train, method='randomForest')
# sc.obj.train <- make.predictions(sc.obj.train)
# sc.obj.train <- evaluate.predictions(sc.obj.train)               
params <- colnames(read_tsv("evaluate_input.txt"))   
params <- map_chr(params, function(x) args[[x]])                                 
#save(sc.obj.train, sc.obj.test, file =  str_c(args[['BASEDIRG']], "results", str_c(args[['MLalgorithm']], args[['lowAbFiltering']], args[['condition']], args[['block']], args[['subsample']], args[['genusLevel']], args[['featureEng']], ".rsave", sep = "___"), sep = "/"))
print(params)                
save(sc.obj.train, sc.obj.test, file = str_c(args[['BASEDIRG']], "results", str_c(str_c(params, collapse = "___"), ".rsave"), sep = "/"))
