This subfolder contain the pretrained REDInet model and the python script to perform inferences on REDItools tabix-indexed output tables. 
It'suggested to operate in an Anaconda or Miniconda virtual environment.
To install Anaconda or Miniconda follows the instructions in: 
    1) https://docs.anaconda.com/free/anaconda/install/linux/
    2) https://docs.anaconda.com/free/miniconda/miniconda-install/
To create a virtual environment run the following commands:
    conda create --name <environment name> python=3.9.0
Samtools must be installed in the virtual environment. To do so run the following commands:
    1) conda activate <environment name>
    2) conda install bioconda::samtools (or conda install bioconda/label/cf201901::samtools) 
    3) conda deactivate
To proper index the tables with tabix run the following command:
    1) conda activate <environment name> 
    2) tabix -s 1 -b 2 -e 2 -c R <output table name>
    3) conda deactivate
REDInet_inference.py script usage require the installation of six python libraries in an python3 environment. To install the required libraries using pip enter the following commands:
    1) conda activate <environment name>
    2) pip install numpy pandas pysam scikit-learn tqdm tensorflow[and-cuda]==2.14.0 
    3) conda deactivate
To perform inference on sequences with gaps it's necessary to download a reference genome il fasta file; only GRCh37 or GRCh38 assembly are allowed.
To download the reference genomes enter the Utilities folder of this package and run these command:
    1) wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
    2) wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
The reference genomes must be decompresed and indexed with samtools using the following commands:
    1) conda activate <environment name>
    1) bgzip -d <compressed reference genome file name>
    1) samtools faidx <reference genome fasta file name>
    3 <reference genome fasta file name> must be: GRCh38.primary_assembly.genome.fa or GRCh37.primary_assembly.genome.fa
To make prediction on the REDItools tabix-indexed output tables with the REDInet.h5 model run the following commands:
    1) conda activate <environment name>
    2) python3 REDInet_Inference.py <OPTIONS>
