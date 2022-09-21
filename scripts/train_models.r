source('/g/scb2/zeller/karcher/dysb_classif/scripts/packages.r')
args <- parseArgs()
load(str_c(args[['BASEDIRG']], "data", "siamcat", "profiles_merged_with_metadata.rimage", sep = "/"))

# Low-prevalence filtering
if (args[['lowAbFiltering']]) {
    print("Doing low abundance filtering...")
    filterBool <- apply(dataMetaAdaptedTruncWideAll, 1, function(x) mean(x>0) > 0.05)
    dataMetaAdaptedTruncWideAll <- dataMetaAdaptedTruncWideAll[filterBool, ]
    stopifnot(all(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS)))   
} else {
    print("NOT doing low abundance filtering...")    
}
                        
if (args[['subsample']]) {
    print("Subsampling datasets....")
    set.seed(1337)
    metaDataWGS <- metaDataWGS %>% 
        mutate(sampleID = rownames(.)) %>%
        group_by(condition) %>% nest() %>% 
        mutate(data = map2(data, condition, function(x, co) {
    # condition here can contain stuff like 'UC/CD' as well as 'UC'. UC are only cases, so the resampling is right now a bit funny.
    if (co == "No condition") {
        return(x)
    }
    #print('Here.')       
    if (dim(x)[1] > 200) {
        return(x %>%     
        sample_n(200, replace = F))
    } else {
        return(x)
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
} else {
    print("NOT subsampling datasets...")
}
                        
if (args[['genusLevel']]) {
    print("Transforming profiles to genus level...")
    tax <- read_tsv('../data/motus2.6_taxonomy_NCBI.tax') %>%
    select(mOTUs_ID, Genus, profiled)
    tax <- tax %>% 
        mutate(profile_corrected = map_chr(profiled, function(x) {
            return(str_c(str_split(x, " ")[[1]][2:length(str_split(x, " ")[[1]])], collapse = " "))
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
    stopifnot(all(colnames(dataMetaAdaptedTruncWideAll) == rownames(metaDataWGS)))
} else {
    print("Not transforming profiles to genus level...")    
}
                        
#if (args[['featureEng']]) {   
### I used to calculcate these indices here.
### Now instead they're being calculated in the load_and_merge_data notebook.    
#     allDat <- getFeatures(dataMetaAdaptedTruncWideAll, metaDataWGS) 
#     morePredictorNames <- allDat[[2]]
#     allDat <- allDat[[1]]    
#     nn <- names(allDat)                           
#     stopifnot(all(map_dbl(allDat, length) == dim(dataMetaAdaptedTruncWideAll)[2]))    

#     tmpData <- matrix(unlist(allDat), byrow = T, nrow = length(allDat)) %>%
#                                      as.data.frame()
#     colnames(tmpData) <- colnames(dataMetaAdaptedTruncWideAll)
#     rownames(tmpData) <- nn   
#     tmpData <- t(tmpData)
#     stopifnot(all(rownames(tmpData) == rownames(metaDataWGS)))
#     metaDataWGS <- metaDataWGS %>% 
#                                      as.data.frame() %>% 
#                                      mutate(sampleID = rownames(.)) %>% 
#                                      left_join(tmpData %>% 
#                                                as.data.frame() %>% 
#                                                mutate(sampleID = rownames(.)), 
#                                                by ='sampleID')                                
#     rownames(metaDataWGS) <- metaDataWGS$sampleID
#     metaDataWGS$sampleID <- NULL       
#}                       
                        
cond <- args[['condition']]

meta.train <- metaDataWGS %>% 
    # Never use Diabetes datasets (for now)
    # TODO: Figure out why not?                     
    filter(condition != "T2D")        

if (args[['includeCohortStudies']]) {
    print("Including cohort studies...")
}  else {
    print("NOT including cohort studies...")
    meta.train <- meta.train %>%
    # Remove 'cohort' cohorts...
    filter(condition != "No condition")
}                     

                        
                        
meta.train <- meta.train %>%                   
    # If we evaluate a mixed IBD cohort, for example, we want to train on (potentially) shared controls. Have a look at Franzosa_NatMed_201
    # More generally, if a cohort contains more than one condition, what do we want to do?
    filter(condition != cond | caseControls == "control") %>%
    as.data.frame()       
                                            
print(str_c("Training on everything but ", cond, "(maybe some of their controls...)"))                        
print(dim(meta.train))
print(dim(meta.train %>% filter(condition == "No condition")))                     
                        

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
                                 
if (args[['block']] == 'condition') {
    print("Blocking by condition at training time...")
    sc.obj.train <- create.data.split(sc.obj.train,
    num.folds = 5, num.resample = as.integer(args[['reSample']]), inseparable = 'condition')
} else {
    print("Not blocking by condition at training time...")
    sc.obj.train <- create.data.split(sc.obj.train,
    num.folds = 5, num.resample = as.integer(args[['reSample']]))
}                        
                    

# potentially: Add features
if (args[['featureEng']]) {
    # Comes from a time where we did feature engineering in this script :)
    #sc.obj.train <- add.meta.pred(sc.obj.train, morePredictorNames)
    print("Adding features (at training time)...")
    sc.obj.train <- add.meta.pred(sc.obj.train, colnames(metaDataWGS)[!colnames(metaDataWGS) %in% c("caseControls", "dataset", "condition")])
} else {
    print("NOT adding features (at training time)...")
}                                 
                        
# train model
if (args[['MLalgorithm']] == "RF") {
  print("Training RF model...")
  sc.obj.train <- train.model(sc.obj.train, method='randomForest')  
} else if (args[['MLalgorithm']] == "lasso") {
  print("Training lasso model...")    
  sc.obj.train <- train.model(sc.obj.train, method='lasso')  
} else if (args[['MLalgorithm']] == "lasso_ll") {
  print("Training lasso (libLinear) model...")    
  sc.obj.train <- train.model(sc.obj.train, method='lasso_ll')  
} else {
  print("ML Algorithm undefined. exciting.")
  stop()
}
print("Done training model...")
          
# apply trained models to held-out dataset
meta.test <- metaDataWGS %>% 
    filter(str_detect(condition, cond)) %>%
    as.data.frame()
                        
print(dim(meta.test))
if (dim(meta.test)[1] == 0) {
    print("You dont have testing data... Exiting")
    adsaads
}                        
                        
                        
                        
save.image('/g/scb2/zeller/karcher/test')                        
                        
                        
                        
                        
sc.obj.test <- siamcat(feat=dataMetaAdaptedTruncWideAll[, colnames(dataMetaAdaptedTruncWideAll) %in% rownames(meta.test)], meta=meta.test,
                        label='caseControls', case='case')

# Explicitly normalize features in test set so that you can then add metadata features to test set                                
sc.obj.test <- normalize.features(sc.obj.test, norm.method = 'log.std',
    norm.param=norm_params(sc.obj.train),feature.type = 'original')      
                                 
# Add features
if (args[['featureEng']]) {
    print("Adding features at testing time...")
    sc.obj.test <- add.meta.pred(sc.obj.test, colnames(metaDataWGS)[!colnames(metaDataWGS) %in% c("caseControls", "dataset", "condition")])
}   
                                                          
# make holdout predictions. DONT NORMALIZE AGAIN
sc.obj.test <- make.predictions(sc.obj.train, 
                                siamcat.holdout = sc.obj.test,
                                normalize.holdout = FALSE)                                 
                                                                  
sc.obj.test <- evaluate.predictions(sc.obj.test)
                        
params <- colnames(read_tsv("evaluate_input.txt"))   
params <- map_chr(params, function(x) args[[x]])                                 
print(str_c("Saving files to ", str_c(args[['BASEDIRG']], "results", str_c(str_c(params, collapse = "___"), str_replace_all(str_replace_all(Sys.time(), " ", "__"), ":", "_"),".rsave"), sep = "/")))                  
save(sc.obj.train, sc.obj.test, file = str_c(args[['BASEDIRG']], "results", str_c(str_c(params, collapse = "___"), str_replace_all(str_replace_all(Sys.time(), " ", "__"), ":", "_"),".rsave"), sep = "/"))              
print("Done saving files :D")                  