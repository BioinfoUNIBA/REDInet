# REDInet
A Temporal Convolutional Network based Python package for the detection of A-to-I RNA editing events in RNA sequencing experiments.

## **Virtual environment preparation and required softwares installation**:
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
   
6) Create REDInet dedicated folder and REDInet package download: <br />

       mkdir REDInet
       cd REDInet
       wget https://github.com/BioinfoUNIBA/REDInet/tree/main/Package

7) Dowload and install REDItools 3. Package and installation guide at: <br />

       https://github.com/BioinfoUNIBA/REDItools3

## **REDInet installation**:
REDInet runs on Tabix indexed REDItools output tables. Both REDInet and REDItools should be installed in a dedicated Anacoda or Miniconda environment. 
To create a the required enviroment: 

The action of every agent <br />
  into the world <br />
starts <br />
  from their physical selves. <br />


3) Create a new virtual environment (it's suggested to create, use and activate a base conda environment with all the required software):

		# create a new conda environment
        conda create --name CtoU python=3.8

		# activate the conda env
		conda activate CtoU

		# install samtools
		conda install -c bioconda samtools

		# install minimap2
		conda install -c bioconda minimap2

		# install f5c (optimised re-implementation of the call-methylation and eventalign modules in Nanopolish)
		conda install -c bioconda f5c

		# create virtual environment inside Conda CtoU env
		python3 -m venv venv
 
