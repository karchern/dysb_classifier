Data="/g/scb2/zeller/karcher/dysb_classif/data/siamcat/profiles_merged_with_metadata.rimage"
BASEDIR="/scratch/karcher/dysbiosis"
BASEDIRG="/g/scb2/zeller/karcher/dysb_classif"
#OUTPUTDIR="/scratch/karcher/dysbiosis/output/"
cp $Data $BASEDIR

conda activate r_test

##############
### PARAMS ###
# This script generates a txt file of input instructions
Rscript get_input_combinations.r 
cat cluster_input.txt | sed "s%$% --BASEDIR=$BASEDIR%" | sed "s%$% --BASEDIRG=$BASEDIRG%"  > tmp.txt
mv tmp.txt cluster_input.txt

echo "Sending the following off to sbatch..."
head cluster_input.txt
sbatch train_models_sbatch.sh cluster_input.txt 
