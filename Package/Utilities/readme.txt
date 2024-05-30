This subfolder contain the pretrained REDInet model and the python script to perform inferences on REDItools tabix-indexed output tables. 
It'suggested to operate in an Anaconda or Miniconda virtual environment.
To install Anaconda or Miniconda follows the instructions in: 
     https://docs.anaconda.com/free/anaconda/install/linux/
     https://docs.anaconda.com/free/miniconda/miniconda-install/
To create a virtual environment run the following commands:
    conda create --name <environment name> python=3.9.0
Samtools must be installed in the virtual environment. To do so run the following commands:
     conda activate <environment name>
     conda install bioconda::samtools (or conda install bioconda/label/cf201901::samtools) 
     conda deactivate
To proper index the tables with tabix run the following command:
     conda activate <environment name> 
     tabix -s 1 -b 2 -e 2 -c R <output table name>
     conda deactivate
REDInet_inference.py script usage require the installation of six python libraries in an python3 environment. To install the required libraries using pip enter the following commands:
     conda activate <environment name>
     pip install numpy pandas pysam scikit-learn tqdm tensorflow[and-cuda]==2.14.0 
     conda deactivate
To perform inference on sequences with gaps it's necessary to download a reference genome il fasta file; only GRCh37 or GRCh38 assembly are allowed.
To download the reference genomes enter the Utilities folder of this package and run these command:
     wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
     wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
The reference genomes must be decompresed and indexed with samtools using the following commands:
     conda activate <environment name>
     bgzip -d <compressed reference genome file name>
     samtools faidx <reference genome fasta file name>
<reference genome fasta file name> must be: GRCh38.primary_assembly.genome.fa or GRCh37.primary_assembly.genome.fa
To know all the REDInet_Inference.py available options run the following command:
     conda activate <environment name>
     python3 REDInet_Inference.py -h
To make prediction on the REDItools tabix-indexed output tables with the REDInet.h5 model run the following commands:
     conda activate <environment name>
     python3 REDInet_Inference.py <OPTIONS>
