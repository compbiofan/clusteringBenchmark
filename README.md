# clusteringBenchmark
Benchmark of the clustering of scDNAseq cells. 

## Table of Contents
- [Commands to run BnpC, SCG and SCClone.](#commands_3methods)
    * [BnpC](#bnpc)
    * [SCG](#scg)
    * [SCClone](#scclone)

# <a name="usage_of_scDNAseq_clustering"></a>Usage of scDNAseq Clustering.
## <a name="software_requirements"></a>Software Requirements ##

1. Python 2.7.15 or up.
2. Python modules: numpy, pandas, pyYaml, scipy.

## <a name="environment_setup"></a>Environment Setup ##

# <a name="commands_3methods"></a>Commands to run BnpC, SCG and SCClone. #

## <a name="bnpc"></a>BnpC ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download BnpC from https://github.com/cbg-ethz/BnpC. There are two ways to use BnpC: 
	* Follow the instructions in their GitHub page and install their conda environment.
	* Use the run_BnpC.py script to run BnpC in the following way:
	``` time python ../../BnpC/run_BnpC.py inputFile -pp 0.75 0.75 -o outputFile ```

	We used ``` pp ``` argument because initially our date had few mutations and we followed the NOTE given by BnpC authors:

	If you run BnpC on panel data with few mutation only or on error free data, we recommend changing the ``` -pp ``` argument to beta distribution closer to uniform, like ``` -pp 0.75 0.75 ``` or even ``` -pp 1 1 ```. Otherwise, BnpC will incorrectly report many singleton clusters.

2. BnpC input file D matrix has to be a binary matrix in the format of rows = mutations and columns = cells. The matrix entries are defined as following: 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating missing value. Any normal D matrix file with dimensions cells x mutations can be 
converted into BnpC file format with the following script:
	``` python processSimInput.py -input inputDFile -output bnpcDFile ```

	OR 
	use ``` -t ``` flag to let BnpC know to transpose the input matrix.

## <a name="scg"></a>SCG ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download SCG from https://github.com/Roth-Lab/scg. Since running SCG requires using multiple user parameters we recommend using our script to run SCG making sure that the path to SCG is valid.
	* Use the following script to run SCG:
	``` python save_multipleSCGresults.py -opDir outputDir -input inputFile -scg_config SCG_CONFIG_FILE_PATH -config_path PATH_TO_SAVE_NEW_CONFIG -sim true/false -niters 10000000 ```
	* Parameters: 
		** ``` -opDir ```, path to save the SCG results.
		** ``` -input ```, input D matrix.
		** ```  -scg_config ```, SCG config file which the script reads and update to use it for running SCG.
		** ``` -config_path ```, path to save the updated config file to run SCG.
		** ``` -sim ```, pass it as True if the input matrix if from Simulated data otherwise pass it as False.
		** ``` -niters ```, number of iterations to use for SCG.
	* Our script restarts SCG for 20 times and based on the MAX_ELBO value choosed the best seed value. Then re-run SCG with this seed value to get the results.

2. SCG input D matrix needs to have cell_IDs annotated and a header indicating mutations with dimensions cells x mutations. This can be achieved by the following script:
	``` python processSimInput.py -input inputDFile -output SCG_DFile ```
	The input matrix should be a binary matrix with 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating missing value. 

## <a name="scclone"></a>SCClone ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download SCClone from https://github.com/qasimyu/scclone. There are two ways to use SCClone:
	* Follow the instructions in their GitHub page and build their binary.
	* Use ``` scclone ``` script inside the scclone folder to run in the following way:
	``` time ../../scclone-1.0/bin/scclone -i inputFile -a 0.01 -b 0.1 -o outputFile ```
	* Parameters:
		** ``` -a ```, false positive rate and default is 0.01.
		** ``` -b ```, false negative rate. If not mentioned then SCClone does a grid search to find the optimum beta value.

2. SCClone accept input D matrix with or withour header. It accepts a binary genotype matrix with cells as rows and mutations as columns. The columns are tab separated with 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating genotype information is missing.


