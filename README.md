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
   
7) Create REDInet dedicated folder and REDInet package download: <br />

       mkdir REDInet
       cd REDInet
       wget https://github.com/BioinfoUNIBA/REDInet/tree/main/Package

## **REDInet Basic Usage**:
The C_to_U_classifier is able to ameliorate the C-to-U editing signal in direct-RNA Nanopore runs. 
It does this by mean of two main pipelines:
