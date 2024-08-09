# REDInet
A Temporal Convolutional Network based Python package, inspired by Google DeepMind' WaveNet, for A-to-I RNA editing events detection in RNA sequencing experiments.

## **Virtual Environment Preparation and Software Installation**:
REDInet runs on Tabix indexed REDItools output tables. <br />
Both REDInet and REDItools should be installed in a dedicated Anacoda or Miniconda environment. <br />
To create the required enviroment:
1) Install Anaconda or Miniconda following the instructions at: <br />
  Anaconda <br />
            
            https://docs.anaconda.com/free/anaconda/install/linux 
   
     Miniconda <br />
            
            https://docs.anaconda.com/free/miniconda/miniconda-install
2) Create the virtual environment with Python version 3.9.0: <br />

       conda create --name REDInet python=3.9.0
   
3) Activate the virtual environment: <br />

       conda activate REDInet

4) Install Samtools: <br />

       conda install bioconda::samtools 

5) Install python required packages: <br />

       pip install numpy==1.26.4
       pip install pandas==2.2.2
       pip install pysam==0.22.1
       pip install scikit-learn==1.4.2
       pip install tqdm==4.66.2
       pip install tensorflow[and-cuda]==2.14.0

7) Dowload and install REDItools. Package and installation guide at: <br />
   REDItools <br />

           https://github.com/BioinfoUNIBA/REDItools
   
   REDItools2 <br />

           https://github.com/BioinfoUNIBA/REDItools2
   
   REDItools3 <br />

           https://github.com/BioinfoUNIBA/REDItools3
   
8) Create a REDInet dedicated folder and REDInet package download: <br />

       mkdir REDInet
       cd REDInet
       wget https://github.com/BioinfoUNIBA/REDInet/tree/main/Package
       cd Package
   
9) Download and prepare GRCh37 and GRCh38 reference genomes: <br />

       cd Utilities
   
       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
   
       gunzip GRCh38.primary_assembly.genome.fa.gz
   
       samtools faidx GRCh38.primary_assembly.genome.fa
   
       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
   
       gunzip GRCh37.primary_assembly.genome.fa.gz
   
       samtools faidx GRCh37.primary_assembly.genome.fa
   
10) Make REDInet_inference.py executable:

        chmod u+x REDInet_Inference.py
    
        cd ..

## **REDInet Usage**:
REDInet classifies A to G substitutions, in RNAseq data, as derived from A-to-I RNA editing or not. <br />
The REDInet pipeline requires BAM files to be prepared via one of the three version of the REDItools package. <br />
In the REDInet/Package/Data folder are provided 2 BAM files as example.  <br />
To run REDInet pipeline using REDItools (version 1) follows these steps:   <br />

1) Launch REDItools REDItoolDnaRNA.py script onto the example BAM file with the following setting: <br />

       python ../REDItools/main/REDItoolDnaRna.py -o ../REDInet/Package/Data -i ../REDInet/Package/Data/SRR12492027.SRR12492028.Aligned.sortedByCoord.out.chr10.bam -f ../REDInet/Package/Utilities/GRCh37.primary_assembly.genome.fa -t 40 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p -s 2 -g 2 -S
   
       python ../REDItools/main/REDItoolDnaRna.py -o ../REDInet/Package/Data -i ../REDInet/Package/Data/SRR12492045.SRR12492046.Aligned.sortedByCoord.out.chr10.bam -f ../REDInet/Package/Utilities/GRCh37.primary_assembly.genome.fa -t 40 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p -s 2 -g 2 -S

3) Compress and Tabix indexing the REDItools output tables: <br /> 

       cd ../REDInet/Package/Data/DnaRna_<REDItools SRR12492045.SRR12492046 numeric ID>
    
       bgzip outTable_<REDItools SRR12492045.SRR12492046 numeric ID>
    
       tabix -s 1 -b 2 -e 2 -c R outTable_<REDItools SRR12492045.SRR12492046 numeric ID>.gz
       
       cd ../REDInet/Package/Data/DnaRna_<REDItools SRR12492045.SRR12492046 numeric ID>
   
       bgzip outTable_<REDItools SRR12492045.SRR12492046 numeric ID>
   
       tabix -s 1 -b 2 -e 2 -c R outTable_<REDItools SRR12492045.SRR12492046 numeric ID>.gz

5) Launch REDInet help message (see REDInet Options section below for further details): <br />

       cd ../REDInet/Package/Utilities
       python3 REDInet_Inference.py -h

6) Example of REDInet basic usage on the REDItools output table using the script enabling missing values imputation (using the reference genome: hg38 or hg19; BETA FEATURE):

       python3 ../REDInet/Package/Utilities/REDInet_Inference.py \
          --I <REDItools_outTable_foldepath> # Input REDItools outTable \
          --O <output_folder> # folder where results will be saved \
          --A <assembly> # Build version: GRCh37 or GRCh38 - REQUIRED \
          --C <MinCov> # e.g. 50 - NOT REQUIRED \
          --M <MinAG> # e.g. 3 - NOT REQUIRED \
          --N <REDItools_outTable_name>

7) Alternatively, a stable and lighter version of the REDInet inference script was released when imputation of missing values is not required or there is a need to use a different genome build:
   
       python3 ../REDInet/Package/Utilities/REDInet_Inference_light_ver.py \
          -r <REDItools_outTable_filepath> # the path for the input REDItools outTable.gz compressed file indexed with tabix (see point 2 above) - REQUIRED \
          -o <output_files_prefix> # the full path for the prefix used to write the 2 output files \
          -c <MinCov> # e.g. 50 - NOT REQUIRED \
          -s <MinAG> # e.g. 3 - NOT REQUIRED \
          -f <MinAGfreq> # e.g 0.01 - NOT REQUIRED

## **REDInet output**:
REDInet analysis is a 3 steps process.
At the end of each step a file is automatically produced and stored in the choiced output directory. 
1) Genomic positions are filtered on the basis of minimum base coverage, minimum A to G substitution rate and minimum number of guanosine in place of adenines.  <br />
   Information regarding genomig positions filtered out in this phase are stored in the file: <br /> <br />
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;../<results folder>/<.gz name>_discarded.txt <br /> <br />
   If no genomic position is filtered out than the <.gz name>_discarded.txt file is empty. <br /> <br />
3) Genomic positions filtered in in the previous phase are futher filtered for missing nucleotides in 101 nucleotides regions centered in them.  <br />
   Extracted sequences having nucleotides with unassigned strands, are discarded and their missing nucleotides are described by the symbol "0 *". <br />
   Extracted sequences having nucleotides with unmatching strands, are discarded and their missing nucleotides are described by the symbol "0 **". <br />
   Extracted sequences having nucleotides with unmatching strands, can be retained and corrected using the appropriate REDInet_Inference.py flag. <br />
   Information regarding genomig positions filtered out in this phase are stored in the file: <br /> <br />
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;../<results folder>/<.gz name>_incompleates.txt <br /> <br />
   If there aren't sequences with missing nucleotides than the <.gz name>_incompleates.txt file is empty. <br /> <br />
3) Genomic positions filtered in in the previous phase are subjected to the classification by the pretrained Deep Neural Network model (REDInet.h5 file).  <br />
   Information regarding the classified genomig positions are stored in the file: <br /> <br />
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;../<results folder>/<.gz name>_predictions.txt <br /> <br />
   If no sequence has been produced, due to the previous filtering steps,  than the <.gz name>_predictions.txt file is empty. <br /> 
          
## **REDInet Options**:
REDInet settings can be tuned to accomodate specific analysis needs.  <br />
This is the list of available parameters that can be set: <br />

    --I                Absolute path to the folder containing the files to be
                       analyzed.By default is the Data folder in this package.
                       
    --O                Absolute path to the folder in which the results are
                       going to be stored.By default is the Results folder in
                       this package.
                       
    --C                Minimum number of RNAseq reads mapping in the positions
                       to be analyzed. Value must be a natural number greater
                       than zero. By defalult the value is 50.
                       
    --F                Minimun percentage of A to G substitutions in the positions
                       to be analyzed. Value must be a floating point number
                       grater than zero. By default the value is 0.01.
                       
    --M                Minimum detected guanosines number in place of adenines
                       in the positions to be analyzed. Value must be a natural
                       number greater than zero. By default the value is 3.
                       
    --N                Compleate list of .gz files to be analyzed. List must be
                       included in box brackets and names must be separated by
                       commas without spaces.
                       Es:[data1.gz,data2.gz,data3.gz]. By default the
                       tool analyze all the .gz files in the input directory.
                       
    --P                Choose whether to work with multiprocessing or
                       not. Possible choices are: yes or no. By default is no.
                       
    --G                Choose to accept sequences with maximum 100 missing
                       nucleotides. Possible choices are: yes or no. By default
                       is no.
                       
    --A                Human genome assembly to use in handling missing
                       nucleotides. Possible choices are: GRCh37 or GRCh38. By
                       default is GRCh37.
                       
    --S                Number of missing nucleotides to be imputed. Value must
                       be an integer number grater than zero and equal or
                       smaller than 100. It is suggested not to use numbers higher than 30.
                       By default the value is 10.

    --R                Correct nucleotides strands in each extracted sequence, 
                       to match corresponding identified site srand. If it's set to no,
                       sequences with nucleotides mapping on different strands are discarded.
                       Possible choices are: yes or no. 
                       By default is no. 
                       
    --U                Declare if the original samples cames from Unstranded RNAseq.
                       Possible choices are: yes or no.
                       By default is no.


## **REDInet ligh version Options and Outputs**:

REDInet_inference_light_ver.py script can be used when imputation of missing values is not required or you need to align your reads against different genome builds than hg38 and hg19.
Here the available options:

      python3 ../REDInet/Package/Utilities/REDInet_Inference_light_ver.py -h                                                                                                                             
      usage: REDInet_Inference_light_ver.py [-h] -r REDITABLE [-m MODEL] [-o OUTPUT_TABLE_PREFIX] [-c COV_THRESHOLD] [-f AGFREQ_THRESHOLD] [-s MINAGSUBS]
      
      REDInet_Inference_light_ver.py
      
      optional arguments:
        -h, --help            show this help message and exit
        -r REDITABLE, --reditable REDITABLE
                              --reditable: a <str> with the fullpath for the input reditable file.
        -m MODEL, --model MODEL
                              --model: a <str> with the fullpath for the model file. [REDIportal model, also the prototype model trained on the kindey dataset can be used at the path REDInet/Notebooks&Scripts/CNN_wavenet_kindney_first_model_prototype/cnn_wavenet_14112023/model_WaveNet_small_log_preprocessing14_11_2023_15_01_48.h5]
        -o OUTPUT_TABLE_PREFIX, --output_table_prefix OUTPUT_TABLE_PREFIX
                              --output_table_prefix: a <str> Indicating the prefix for the output files. [None]
        -c COV_THRESHOLD, --cov_threshold COV_THRESHOLD
                              --cov_threshold: a <int> Indicating the minimum coverage to make inference. [50]
        -f AGFREQ_THRESHOLD, --AGfreq_threshold AGFREQ_THRESHOLD
                              --AGfreq_threshold: a <float> Indicating the minimum AG substitution frequency to make inference. [0.01]
        -s MINAGSUBS, --minAGsubs MINAGSUBS
                              --minAGsubs: a <int> Indicating the minimum AG substitutions to make inference. [3]

The script will iterate over the compressed and tabix indexed REDItools outTable and it will produce two different files with the prefix used into the -o option:

A) <output_files_prefix>.feature_vectors.tsv: Tabular file containing features vectors of sites with complete intervals (no missing values) satisfying the selected filters.

B) <output_files_prefix>.predictions.tsv: Tabular files with the predictions made by REDInet selected model (by default it uses the REDIportal model). Each row contains the prediction at a per-site level and has 11 columns:

   1) region ------> the genomic region
   2) position ----> the genomic position
   3) Strand ------> the transcript strand infered by REDItools
   4) FreqAGrna ---> the AG substitution frequency
   5) start -------> the start position of the 101 nt-long interval used to perform the REDInet prediction
   6) stop --------> the stop position of the 101 nt-long interval used to perform the REDInet prediction
   7) int_len -----> the expected length of the interval
   8) TabixLen ----> the actual length of the extracted interval from the tabix indexed REDItools outTable
   9) snp_proba ---> Probability for the current site being a SNP (negative class)
   10) ed_proba ---> Probability for the current site being an Editing Site (positive class)
   11) y_hat  -----> Output class computed via softmax functions on SNP/Editing probabilities. 0: Predicted SNP / 1: Predicted Editing Site


## **Notes**:
REDInet pipeline is optimizes to run on REDItools protocol-derived BAM files. <br />
So it's recommended to produce the BAM files via the REDItools protocol at:  <br />

     https://www.nature.com/articles/s41596-019-0279-7
     
REDInet is set-up to work on stranded RNAseq-derived BAM files. <br />
In case of unstranded RNAseq-derived BAM files, it's required to infer the strand using the choiced assembly GENCODE annotation GTF file.
The GTF file has to be sorted and tabix-indexed before its usage. <br />
This procedure can only be performed with REDItools (version 1). <br /> 
It's suggested to use the REDItoolDnaRnav13.py in the NPscripts REDItools folder. <br />
In this case REDInet_Inference.py should be used using the appropriate flag: --U yes. <br />
REDInet is compatible with every versions of REDItools.  <br />

