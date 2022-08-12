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

The output of BnpC has the following two important files:

1. The cells assigned to each cluster number is provided in the file 'assignment.txt' and the genotype of each of these clusters are present in the file 'genotypes_posterior_mean.tsv'.
2. Using these two files we can get the consensus genotype matrix. Use the following script to get the consensus genotype matrix:

	``` python bnpc_getGmatrix.py -cc assignment.txt -gp genotypes_posterior_mean.tsv -sim true -op consensus_genotype_fileName.tsv ```
3. To evaluate the accuracy, sensitivity, specificity use the following script:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -sim true > eval_metrics.txt ``` 
4. To evaluate the V-measure use the following script:

	``` python ../evaluation.py -i "bnpc:"assignment.txt -G groundTruthFile -v >> eval_metrics.txt ```

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
	* Our script restarts SCG for 20 times and based on the MAX_ELBO value chose the best seed value. Then re-run SCG with this seed value to get the results.

2. SCG input D matrix needs to have cell_IDs annotated and a header indicating mutations with dimensions cells x mutations. This can be achieved by the following script:
	``` python processSimInput.py -input inputDFile -output SCG_DFile ```
	The input matrix should be a binary matrix with 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating missing value. 

The output of SCG has the following two important files:

1. 'cluster_posteriors.tsv.gz' has the probability of clusters belonging to each cell. We have to choose the cluster number with the maximum probability for a cell. 
From the file 'genotype_posteriors.tsv.gz'  we get the genotype for each cluster.
2. Using these two files we can get the consensus genotype matrix. Use the following script to get the consensus genotype matrix:

	``` python scg_getGmatrix_new.py -cp cluster_posteriors.tsv -gp genotype_posteriors.tsv.gz -D inputDMatrix -output consensus_genotype_fileName.tsv ```
3. To evaluate the accuracy, sensitivity, specificity use the following script:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -sim true > eval_metrics.txt ```
4. To evaluate the V-measure use the following script:

        ``` python ../evaluation.py -i "scg:"cluster_posteriors.tsv -G groundTruthFile -v >> eval_metrics.txt ```

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


The output of SCClone has the following two important files:

1. The file 'data.cell_assignment' provides the assigned cluster number for each cell. And the file 'data.clone_genotypes' provided the genotypes for each cluster.
2. Using these two files we can get the consensus genotype matrix. Use the following script to get the consensus genotype matrix:

	``` python scclone_getGmatrix.py -ca data.cell_assignment -cg data.clone_genotypes -sim true -op consensus_genotype_fileName.tsv ```
3. To evaluate the accuracy, sensitivity, specificity use the following script:

        ``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -sim true > eval_metrics.txt ```
4. To evaluate the V-measure use the following script:

        ``` python ../evaluation.py -i "scclone:"data.cell_assignment -G groundTruthFile -v >> eval_metrics.txt ```

## <a name="simulator"></a>Simulator ##

### Steps to run the simulator using different variables. ###

#### make tree structure for different number of leaves. Replace the values of i with desired number of leaves, folders in t*
``` for i in 4 8 16 32; do for j in `seq 1 5`; do python gen_tree.py -F $i -B 0.2 -o t$i/rep$j/tree_${i}_p2.csv; done; done ```

Here the parameter ``` B ``` indicates beta splitting variable whose default value is 0.2.

#### make data for all different number of leaves.  Replace the values of ```i``` with desired number of leaves, folders in t*
``` for i in 4 8 16 32; do for j in `seq 1 5`; do python sim_par.py -f t$i/rep$j/tree_${i}_p2.csv -P t$i/rep$j/input_t${i}_rep${j}; done; done ```

#### make data when alpha (False positive rate) varies. Replace the value of ```a``` with desired false positive rate, folders in a*
``` a=0.001; i=ap001; for j in `seq 1 5`; do python sim_par.py -a $a -f $i/rep$j/tree_8_p2.csv -P $i/rep$j/input_${i}_rep${j}; done; ```

#### make data when beta (False negative rate) varies. Replace the value of ```b``` with desired false negative rate, folders in b*
``` b=0.1; i=bp1; for j in `seq 1 5`; do python sim_par.py -b $b -f $i/rep$j/tree_8_p2.csv -P $i/rep$j/input_${i}_rep${j}; done; ```

#### make data when missing data rate varies. Replace the value of ```m``` with desired missing rate, folders in m*
``` m=0.3; i=mp3; for j in `seq 1 5`; do python sim_par.py -m $m -f $i/rep$j/tree_8_p2.csv -P $i/rep$j/input_${i}_rep${j}; done; ```

#### make data when number of cells varies. Replace the value of ```N``` with desired number of cells, folders in N*
``` N=100; i=N100; for j in `seq 1 5`; do python sim_par.py -c $N -f $i/rep$j/tree_8_p2.csv -P $i/rep$j/input_${i}_rep${j}; done; ```

#### make data when number of mutations varies. Replace the value of ```M``` with desired number of mutations, folders in M*
``` M=50; i=M50; for j in `seq 1 5`; do python sim_par.py -n $M -f $i/rep$j/tree_8_p2.csv -P $i/rep$j/input_${i}_rep${j}; done; ```

#### make tree structure for different beta splitting model variable. Replace the value of ```s``` with desired beta splitting variable, folders in s*
``` s=0.05; i=sp05; for j in `seq 1 5`; do python gen_tree.py -B $s -o $i/rep$j/tree_8_p05.csv; python sim_par.py -f $i/rep$j/tree_8_p05.csv -P $i/rep$j/input_${i}_rep${j}; done ```

#### make data when doublet rate varies. Replace the value of ```dp``` with doublet rate, folders in s*
``` dp=0.1; i=dp1; for j in `seq 1 5`; do python sim_par.py -e $dp -f $i/rep$j/tree_8_p2.csv -P $i/rep$j/input_${i}_rep${j}; done; ```



