library(gridExtra)
library(tidyverse)
# CDI has only controls and thus always fails in evaulation step.
#conditions <- c("CRC","IBD","ACVD","MS","OB","HTN","LCIR","CDI","AS","RA","BC")
#conditions <- c("CRC","IBD","ACVD","MS","OB","HTN","LCIR","AS","RA","BC")
#conditions <- c("CRC","UC","CD","ACVD","MS","OB","HTN","LCIR","AS","RA","BC")
conditions <- c("CRC","ADA", "UC","CD","ACVD","MS","OB","HTN","LCIR","AS","RA","BC")
lowAbFiltering <- c("TRUE", "FALSE")
block <- c("condition", "no_blocking")
MLalgorithm <- c("RF", "lasso", "lasso_ll")
subsample <- c("TRUE", "FALSE")
genusLevel <- c("TRUE", "FALSE")
featureEng <- c("TRUE", "FALSE")

reSample <- c("1", "10")
includeCohortStudies <- c("TRUE", "FALSE")
# batchCorrection <- c("TRUE", "FALSE")

data <- expand.grid(conditions, lowAbFiltering, block, MLalgorithm,subsample, genusLevel, featureEng, reSample, includeCohortStudies, stringsAsFactors = FALSE) %>% 
as.data.frame()
colnames(data) <- c("condition", "lowAbFiltering", 'block', "MLalgorithm", "subsample", "genusLevel", 'featureEng', 'reSample', 'includeCohortStudies')

a <- list()
for (i in 1:dim(data)[1]){
    a[length(a)+1] <-  str_c("Rscript train_models.r ", str_c(str_c("--", colnames(data)), data[i,,drop=T], sep = "=", collapse = " "))
}
data$args <- unlist(a)
data %>% 
select(args) %>%
write_csv("cluster_input.txt", col_names = F)
data %>%
select(-args) %>%
write_tsv("evaluate_input.txt")

