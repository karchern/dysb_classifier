source('/g/scb2/zeller/karcher/dysb_classif/scripts/packages.r')
baseDir <- "/g/scb2/zeller/karcher/dysb_classif/results"
# TODO: expose
args <- parseArgs()
#args['inputFile'] <- "RA___TRUE___no___lasso_ll___TRUE___TRUE___FALSE.rsave"
nrFile <- load(str_c(baseDir, args[['inputFile']], sep = "/"))
print(args)
#######
######
#####
#### This is the correct line. Uncomment me. When you're done.
cns <- colnames(read_tsv("evaluate_input.txt"))
#cns <- c('MLalgorithm', 'lowAbFiltering', 'condition', 'block', 'subsample', 'genusLevel', 'featureEng')
params <- str_split(args[['inputFile']], "/")[[1]][length(str_split(args[['inputFile']], "/")[[1]])]
params <- str_split(params, "___")[[1]]
params <- map_chr(params, function(x) str_replace(x, ".rsave", ""))
#print(params)
#print(cns)                  
stopifnot(length(params) == length(cns))
params <- matrix(params) %>% 
                  t() %>%
                  as.data.frame()
colnames(params) <- cns
params$auc <- sc.obj.test@eval_data$auroc                  
params %>%
                  write_tsv(args[['outputFile']])