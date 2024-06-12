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

       pip install numpy pandas pysam scikit-learn tqdm tensorflow[and-cuda]==2.14.0

6) Dowload and install REDItools 3. Package and installation guide at: <br />

       https://github.com/BioinfoUNIBA/REDItools3
   
7) Create a REDInet dedicated folder and REDInet package download: <br />

       mkdir REDInet
       cd REDInet
       wget https://github.com/BioinfoUNIBA/REDInet/tree/main/Package
       cd Package
   
8) Download and prepare GRCh37 and GRCh38 reference genomes: <br />

       cd Utilities
       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
       bgzib -d GRCh38.primary_assembly.genome.fa.gz
       samtools faidx GRCh38.primary_assembly.genome.fa
       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
       bgzib -d GRCh37.primary_assembly.genome.fa.gz
       samtools faidx GRCh37.primary_assembly.genome.fa
   
9) Make REDInet_inference.py executable:

        chmod u+x REDInet_Inference.py
        cd ..

## **REDInet Usage**:
REDInet classifies A to G substitutions, in RNAseq data, as derived from A-to-I RNA editing or not. <br />
The REDInet pipeline requires BAM files to be prepared via REDItoolDnaRNA.py script from the REDItools3 package. <br />
In the ../REDInet/Package/Data folder are provided 2 BAM files, as example to run REDInet pipeline as follows:  <br />
1) Launch the REDItoolDnaRNA.py onto the example BAM file with the following setting: <br />

       cd Data
       python3 ../REDItools3/main/REDItoolDnaRna.py -o ../REDInet/Package/Data/outTable_SRR12492027.SRR12492028 -i ../REDInet/Package/Data/SRR12492027.SRR12492028.Aligned.sortedByCoord.out.chr10.bam -f ../REDInet/Package/Utilities/GRCh37.primary_assembly.genome.fa -t 40 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p -s 2 -g 2 -S
       python3 ../REDItools3/main/REDItoolDnaRna.py -o ../REDInet/Package/Data/outTable_SRR12492045.SRR12492046 -i ../REDInet/Package/Data/SRR12492045.SRR12492046.Aligned.sortedByCoord.out.chr10.bam -f ../REDInet/Package/Utilities/GRCh37.primary_assembly.genome.fa -t 40 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p -s 2 -g 2 -S

2) Compress and Tabix indexing the REDItools3 output tables: <br /> 
            
       bgzip outTable_SRR12492027.SRR12492028
       tabix -s 1 -b 2 -e 2 -c R outTable_SRR12492027.SRR12492028.gz
       bgzip outTable_SRR12492045.SRR12492046
       tabix -s 1 -b 2 -e 2 -c R outTable_SRR12492045.SRR12492046.gz
       cd ..
3) Launch REDInet analysis on the REDItools3 output table: <br />

       cd Utilities
       python3 REDInet_Inference.py  

## **REDInet output**:
REDInet analysis is a 3 steps process.
At the end of each step a file is automatically produced and stored in the choiced output directory. 
1) Genomic positions are filtered on the basis of minimum base coverage, minimum A to G substitution rate and minimum number of guanosine in place of adenines.  <br />
   Information regarding genomig positions filtered out in this phase are stored in the file: <br /> <br />
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;../<results folder>/<.gz name>_discarded.txt <br /> <br />
3) Genomic positions filtered in in the previous phase are futher filtered for missing nucleotides in 101 nucleotides regions centered in them.  <br />
   Information regarding genomig positions filtered out in this phase are stored in the file: <br /> <br />
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;../<results folder>/<.gz name>_incompleates.txt <br /> <br />
3) Genomic positions filtered in in the previous phase are subjected to the classification by the pretrained Deep Neural Network Model (REDInet.h5 file).  <br />
   Information regarding the classified genomig positions are stored in the file: <br /> <br />
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;../<results folder>/<.gz name>_predictions.txt <br /> 
          
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
                       
    --G                Choose to accept sequences with maximum 30 missing
                       nucleotides. Possible choices are: yes or no. By default
                       is no.
                       
    --A                Human genome assembly to use in handling missing
                       nucleotides. Possible choices are: GRCh37 or GRCh38. By
                       default is GRCh37.
                       
    --S                Number of missing nucleotides to be imputed. Value must
                       be an integer number grater than zero and equal or
                       smaller than 30. By default the value is 10.

## **Notes**:
REDInet pipeline is optimizes to run on REDItools protocol-derived BAM files. <br />
So it's recommended to produce the BAM files via the REDItools protocol at:  <br />

     https://www.nature.com/articles/s41596-019-0279-7

REDInet is also compatible with older versions of REDItools.  <br />

