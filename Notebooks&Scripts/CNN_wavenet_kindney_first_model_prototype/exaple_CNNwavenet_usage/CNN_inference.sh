#!/bin/bash
### Author: Adriano Fonzino email: adriano.fonzino@uniba.it

source /lustrehome/afonzino/anaconda3/bin/activate CtoU
source /lustrehome/afonzino/NanoSpeech/venv/bin/activate

# input
cellline=$1
sample=$2
genotype=$3
assembly=$4
outtable_foldepath=$5
outtable_name=$6
MinCov=$7
MinAG=$8
reditable_path=$outtable_foldepath/$outtable_name

# some logging
echo INPUT:
echo Cell-line: $cellline
echo Sample: $sample
echo Genotype: $genotype
echo Assembly: $assembly
echo REDItools Output Folder: $outtable_foldepath
echo Compressed REDItools OutTable name: $outtable_name
echo Minumum Coverage: $MinCov
echo Compressed REDItools OutTable filepath: $reditable_path

# output paths
outputprefix=/lustre/bio_running/new_basecaller/REDINET_TEST_30_07_2024/REDInet/Package/Results_CNN_inference/$cellline/$sample.$outtable_name

echo OUTPUT:
echo Output filepath: $outputprefix

echo Launching CNN_inference.py with log-preprocessing:
# modify path to h5 model and python executable
python3 /lustre/bio_running/A_to_I_Pietro/model_test_2023/src/kidney_dataset/CNN_inference.py \
    -r $reditable_path \
    -m /lustre/bio_running/A_to_I_Pietro/model_test_2023/src/kidney_dataset/cnn_wavenet_14112023/model_WaveNet_small_log_preprocessing14_11_2023_15_01_48.h5 \
    -o $outputprefix \
    -c $MinCov \
    -f 0.01