# REDInet
A Temporal Convolutional Network based Python package, inspired by Google DeepMind' WaveNet, A-to-I RNA editing events detection in RNA sequencing experiments.

## **Virtual Environment Preparation and Software Installation**:
REDInet runs on Tabix indexed REDItools output tables. Both REDInet and REDItools should be installed in a dedicated Anacoda or Miniconda environment. 
To create a the required enviroment:
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

       conda install bioconda::samtools (or conda install bioconda/label/cf201901::samtools)

5) Install python required packages: <br />

       pip install numpy pandas pysam scikit-learn tqdm tensorflow[and-cuda]==2.14.0

6) Dowload and install REDItools 3. Package and installation guide at: <br />

       https://github.com/BioinfoUNIBA/REDItools3
   
7) Create a REDInet dedicated folder and REDInet package download: <br />

       mkdir REDInet
       cd REDInet
       wget https://github.com/BioinfoUNIBA/REDInet/tree/main/Package
       cd Package
   
9) Download and prepare GRCh37 and GRCh38 reference genomes: <br />

       cd Utilities
       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
       bgzib -d GRCh38.primary_assembly.genome.fa.gz
       samtools faidx GRCh38.primary_assembly.genome.fa
       wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
       bgzib -d GRCh37.primary_assembly.genome.fa.gz
       samtools faidx GRCh37.primary_assembly.genome.fa
       cd ..
       
## **REDInet Usage**:
REDInet classifies A to G substitutions, in RNAseq data, as derived from A-to-I RNA editing or not.  
The basic REDInet pipeline requires BAM files to be prepared via REDItoolDnaRNA.py script from the REDItools3 package. 
In the ../REDInet/Package/Data folder is provided an example BAM file to run REDInet basic pipeline as follows. 
1) Launch the REDItoolDnaRNA.py onto the example BAM file with the following setting: <br />

       cd Data
       python3 ../REDItools3/main/REDItoolDnaRna.py -o ../REDInet/Package/Data/outTable_SRR12492027.SRR12492028 -i ../REDInet/Package/Data/SRR12492027.SRR12492028.Aligned.sortedByCoord.out.chr10.bam -f ../REDInet/Package/Utilities/GRCh37.primary_assembly.genome.fa -t 40 -c 0,1 -m 0,255 -v 1 -q 0,30 -e -n 0.0 -N 0.0 -u -l -p -s 2 -g 2 -S

2) Compress and Tabix indexing the REDItools3 output table: <br /> 
            
       bgzip outTable_SRR12492027.SRR12492028
       tabix -s 1 -b 2 -e 2 -c R outTable_SRR12492027.SRR12492028.gz
            
    

