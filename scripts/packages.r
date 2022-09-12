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