# REDInet
A Temporal Convolutional Network based Python package for the detection of A-to-I RNA editing events in RNA sequencing experiments.

## **Required Softwares**:
REDInet runs on REDItools output tables. Both REDInet and REDItools should be installed in a dedicated Anacoda or Miniconda environment. 
To create a the required enviroment:
1) Install Anaconda or Miniconda following the instructions: <br />
  Anaconda <br />
            
            https://docs.anaconda.com/free/anaconda/install/linux  <br />
   
  Miniconda <br />
            
            https://docs.anaconda.com/free/miniconda/miniconda-install 
  

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
 
