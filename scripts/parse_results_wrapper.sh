#!/bin/bash -e
for file in $(ls ../results | sed "s/\.rsave//")
do 
    echo -e "Rscript parse_results.r --inputFile=../results/${file}.rsave --outputFile=../results_parsed/${file}.tsv"
done | parallel -j 12 {}
