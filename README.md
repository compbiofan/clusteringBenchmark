# clusteringBenchmark
Benchmark of the clustering of scDNAseq cells. 

## Table of Contents
- [Commands to run BnpC, SCG, SCClone, RobustClone, SCITE, SBMClone and simulators.](#commands_3methods)
    * [BnpC](#bnpc)
    * [SCG](#scg)
    * [SCClone](#scclone)
    * [RobustClone](#rcclone)
    * [SCITE](#scite)
    * [SBMClone](#sbmclone)
    * [UltraLow coverage simulator](#ucsim)
    * [Simulator](#simulator)

# <a name="usage_of_scDNAseq_clustering"></a>Usage of scDNAseq Clustering.
## <a name="software_requirements"></a>Software Requirements ##

1. Python 2.7.15 or up.
2. Python modules: numpy, pandas, pyYaml, scipy.

## <a name="environment_setup"></a>Environment Setup ##

# <a name="commands_3methods"></a>Commands to run BnpC, SCG, SCClone, RobustClone, SCITE, SBMClone and simulators. #

## <a name="bnpc"></a>BnpC ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download BnpC from https://github.com/cbg-ethz/BnpC. There are two ways to use BnpC: 
	* Follow the instructions in their GitHub page and install their conda environment.
	* Use the run_BnpC.py script to run BnpC in the following way:
	``` time python BnpC/run_BnpC.py inputFile -o outputFile ```

	We used ``` pp ``` argument because initially our date had few mutations and we followed the NOTE given by BnpC authors:

	If you run BnpC on panel data with few mutation only or on error free data, we recommend changing the ``` -pp ``` argument to beta distribution closer to uniform, like ``` -pp 0.75 0.75 ``` or even ``` -pp 1 1 ```. Otherwise, BnpC will incorrectly report many singleton clusters.

2. BnpC input file D matrix has to be a binary matrix in the format of rows = mutations and columns = cells. The matrix entries are defined as following: 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating missing value. Any normal D matrix file with dimensions cells x mutations can be 
converted into BnpC file format with the following script:
	``` python processSimInput.py -input inputDFile -output bnpcDFile ```

	OR 
	use ``` -t ``` flag to let BnpC know to transpose the input matrix.

We have included an example input data to use which is available at ``` bnpc/example_data/input_D.csv```.

The output of BnpC has the following two important files:

1. The cells assigned to each cluster number is provided in the file 'assignment.txt' and the genotype of each of these clusters are present in the file 'genotypes_posterior_mean.tsv'.
2. Using these two files we can get the consensus genotype matrix. Use the following script to get the consensus genotype matrix:

	``` python bnpc_getGmatrix.py -cc assignment.txt -gp genotypes_posterior_mean.tsv -op consensus_genotype_fileName.tsv ```
3. To evaluate the accuracy, sensitivity, specificity use the following script:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false > eval_metrics.txt ```
The flag ```header``` is used to indicate if the consensus_genotype_fileName.tsv has a header.

	For doublets run the script with a ```doublet``` flag:
	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false -doublet true -doubletFile doubletFileInformation > eval_metrics.txt ```
4. To evaluate the V-measure use the following script:

	``` python ../evaluation.py -i "bnpc:"assignment.txt -G groundTruthFile -v >> eval_metrics.txt ```
For doublets run the script with a ```doublet``` flag:

	``` python ../evaluation.py -i "bnpc:"assignment.txt -G groundTruthFile -d true -df doubletFileInformation -v >> eval_metrics.txt ```
	
	An example doublet information file is provided in ```bnpc/example_data```

## <a name="scg"></a>SCG ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download SCG from https://github.com/Roth-Lab/scg. Since running SCG requires using multiple user parameters we recommend using our script to run SCG making sure that the path to SCG is valid.
	* Use the following script to run SCG:
	``` python save_multipleSCGresults.py -opDir outputDir -input inputFile -scg_config SCG_CONFIG_FILE_PATH -config_path PATH_TO_SAVE_NEW_CONFIG -config_fname NEW_CONFIG_FILENAME -niters 10000000 ```
	* Parameters: 
		* ``` -opDir ```, path to save the SCG results. 
		* ``` -input ```, input D matrix.
		* ```  -scg_config ```, original SCG config file which the script reads and update to use it for running SCG.
		* ``` -config_path ```, path to save the updated config file to run SCG.
		* ``` -config_fname ```, filename that will be used to save the new updated config file to run SCG.
		* ``` -niters ```, number of iterations to use for SCG.
	* Our script restarts SCG for 20 times and based on the MAX_ELBO value chose the best seed value. Then re-run SCG with this seed value to get the results.
We included an example config file in 'scg/example_data' where the parameters ``` num_clusters ```, ``` gamma_prior ```, ``` state_prior ``` and ``` kappa_prior ``` are mentioned. Please refer to their GitHub page for detailed description of all the parameters.

2. SCG input D matrix needs to have cell_IDs annotated and a header indicating mutations with dimensions cells x mutations. This can be achieved by the following script:
	``` python processSimInput.py -input inputDFile -output SCG_DFile ```
	The input matrix should be a binary matrix with 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating missing value. 

We have included a processed example input data to use which is available at ``` scg/example_data/input_D.tsv.gz```.

The output of SCG has the following two important files:

1. 'cluster_posteriors.tsv.gz' has the probability of clusters belonging to each cell. We have to choose the cluster number with the maximum probability for a cell. 
From the file 'genotype_posteriors.tsv.gz'  we get the genotype for each cluster.
2. Using these two files we can get the consensus genotype matrix. Use the following script to get the consensus genotype matrix:

	``` python scg_getGmatrix.py -cp cluster_posteriors.tsv -gp genotype_posteriors.tsv.gz -D inputDMatrix -output consensus_genotype_fileName.tsv ```
3. To evaluate the sensitivity, specificity use the following script:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false > eval_metrics.txt ```
The flag ```header``` is used to indicate if the consensus_genotype_fileName.tsv has a header.

For doublets run the script with a ```doublet``` flag:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false -doublet true -doubletFile doubletFileInformation > eval_metrics.txt ```
	
4. To evaluate the V-measure use the following script:

        ``` python ../evaluation.py -i "scg:"cluster_posteriors.tsv -G groundTruthFile -v >> eval_metrics.txt ```
 For doublets run the script with a ```doublet``` flag:
 
	``` 
	python ../evaluation.py -i "scg:"cluster_posteriors.tsv -G groundTruthFile -d true -df doubletFileInformation -v >> eval_metrics.txt 
	```
An example doublet information file is provided in ```scg/example_data```

## <a name="scclone"></a>SCClone ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download SCClone from https://github.com/qasimyu/scclone. There are two ways to use SCClone:
	* Follow the instructions in their GitHub page and build their binary.
	* Use ``` scclone ``` script inside the scclone folder to run in the following way:
	``` time scclone-1.0/bin/scclone -i inputFile -a 0.01 -o outputFile ```
	* Parameters:
		* ``` -a ```, false positive rate and default is 0.01.
		* ``` -b ```, false negative rate. If not mentioned then SCClone does a grid search to find the optimum beta value.

2. SCClone accept input D matrix with or withour header. It accepts a binary genotype matrix with cells as rows and mutations as columns. The columns are tab separated with 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating genotype information is missing.

We have included an example input data to use which is available at ``` scclone/example_data/input_D.csv```.

The output of SCClone has the following two important files:

1. The file 'data.cell_assignment' provides the assigned cluster number for each cell. And the file 'data.clone_genotypes' provided the genotypes for each cluster.
2. Using these two files we can get the consensus genotype matrix. Use the following script to get the consensus genotype matrix:

	``` python scclone_getGmatrix.py -ca data.cell_assignment -cg data.clone_genotypes -op consensus_genotype_fileName.tsv ```
3. To evaluate the sensitivity, specificity use the following script:

        ``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false > eval_metrics.txt ```
The flag ```header``` is used to indicate if the consensus_genotype_fileName.tsv has a header.

For doublets run the script with a ```doublet``` flag:

	``` 
	python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false -doublet true -doubletFile doubletFileInformation > eval_metrics.txt 
	```
	
4. To evaluate the V-measure use the following script:

        ``` python ../evaluation.py -i "scclone:"data.cell_assignment -G groundTruthFile -v >> eval_metrics.txt ```
For doublets run the script with a ``` doublet ``` flag:
	``` python ../evaluation.py -i "scclone:"data.cell_assignment -G groundTruthFile -d true -df doubletFileInformation -v >> eval_metrics.txt ```
	An example doublet information file is provided in ```scclone/example_data```.

## <a name="rc"></a>RobustClone ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. The folder ``` RobustClone ``` have all the modified scripts required to run the experiments. It accepts a binary genotype matrix cells as rows and mutations as columns. The columns are comma separated with 0 indicating absence of mutation, 1 indicating presence of mutation and 3 indicating genotype information is missing. RobustClone requires its columns and indices to be annotated. We used this script to annotate it ``` python processInput.py -input input.D.csv -output input_rc.csv -mut 200 ``` where we passed the input file, desired output file and the number of mutations as parameters.

We have included an example input data to use available at ``` RobustClone/exampledata/input_rc.csv ```.

2. RobustClone can be run using the MATLAB script 'carryout_RPCA.m' and R script 'carryout_clonal_tree_1.R'. We pass the output file names 'clone_cells.csv' and 'clone_genotype.csv' to the R script. 

	* The scripts should be run in the following order:
	``` 
	matlab -nodisplay -nosplash -nodesktop -r \"try, carryout_RPCA input_rc.csv output_file.mat, catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit\ 
	Rscript carryout_clonal_tree_1.R output_file.mat clone_cells.csv clone_genotype.csv
	```
	* Parameters:
                * Replace 'input_rc.csv' with the input genotype matrix. You can use the example data provided to test.
                * Replace 'output_file.mat' with the output mat file name where you want to save the RPCA results. 
                * Replace 'clone_cells.csv' with the file name where you want to save the cells belonging to cluster result.
                * Replace 'clone_genotype.csv' with the file name where you want to save the clonal genotype.

The output of RobustClone has following two important files which we use to get the consensus genotype and cells assigned to each cluster for evaluation:

1. The file 'clone_cells.csv' where each row indicates a clone and the comma separated row values indicates the cells present in the clones. The file 'clone_genotype.csv' gives us the clonal genotype.
2. We get the consensus genotype and cells assigned to each cluster('assignment.txt') for our evaluation using the following script:

	``` python processOutput.py -cells clone_cells.csv -cloneg clone_genotype.csv -noOfCells 500 -op_ass assignment.txt -output consensus_genotype_fileName.tsv ```
3. To evaluate the sensitivity, specificity use the following script:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true/false > eval_metrics.txt ```
The flag ```header``` is used to indicate if the consensus_genotype_fileName.tsv has a header.

For doublets run the script with a ```doublet``` flag:

	``` python ../evaluateMetrics.py -cg consensus_genotype_fileName.tsv -gtG groundTruthFile -header true -doublet true -doubletFile doubletFileInformation > eval_metrics.txt ```
4. To evaluate the V-measure use the following script:

	``` python ../evaluation.py -i "robustclone:"assignment.txt -G groundTruthFile -v ```
For doublets run the script with a ``` doublet ``` flag:
	``` python ../evaluation.py -i "robustclone:"assignment.txt -G groundTruthFile -d true -df doubletFileInformation -v ```
	An example doublet information file is provided in ```RobustClone/exampledata```.

## <a name="rc"></a>SCITE ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download SCITE from https://github.com/cbg-ethz/SCITE. To compile the C/C++ program, open a terminal and go to the folder containing the source files, and type the following for faster execution of SCITE
	``` clang++ *.cpp -o scite -O3 ```
2. SCITE's input requires the matrix to have the mutations as rows and cells as columns. The genotype values are separated by a space where 1 indicates presence of mutation, 0 indicates absence of mutation and 3 indicates missing data. We used the following script to get such input D matrix ``` python processInput.py -input input.D.csv -output outfileFileName ```. 
3. Once you have SCITE installed you can run it using the following command:
	* ``` SCITE/scite -i example.D.csv -n 200 -m 500 -r 1 -l 1565000 -fd 0.01 -ad 0.1 -a -max_treelist_size 1 -o outputDir ``` 
	* Parameters:
		* ```-i``` indicates the input D matrix to pass. SCITE's input matrix should have mutations as rows and cells as columns. The genotype values are separated by a space where 1 indicates presence of mutation, 0 indicates absence of mutation and 3 indicates missing data.
		* ```-n``` indicates number of mutations. (REQUIRED)
		* ```-m``` indicates number of cells. (REQUIRED)
		* ```-r``` is the desired number of repetitions of the MCMC. (REQUIRED)
		* ```-l``` is the desired chain length of each MCMC repetition. We calculated this value using the formula 10n^2(log n). (REQUIRED)
 		* ```-fd``` is the estimated false positive rate. (REQUIRED)
		* ```-ad``` is the estimated false negative rate. (REQUIRED)
		* ```-a``` with this option SCITE adds the individual cells as additional nodes (leafs) to the reported trees. (OPTIONAL)
		* ```-max_treelist_size``` This limits the number of co-optimal trees written to output files to the specified value. (OPTIONAL)
		* ```-o``` indicates the output directory to save the results. (OPTIONAL)

4. Running SCITE will give an output file having the tree and cells attached to the tree saved as a '.gv' file. From here we can get the clones having the cells and the consensus genotype using the following script:
	* ``` python processOutput.py -input example.D.csv -tree output_ml0.gv -noOfCells 500 -op_ass assignment.txt -op_gen consensus_genotype.csv ``` 
	* Parameters (ALL REQUIRED):
		* ```-input``` pass the input D matrix.
		* ```-tree``` pass the output tree in .gv file.
		* ```-noOfCells``` pass the number of cells in the D matrix.
		* ```-op_ass``` the filename to save the assignment of cells to clones.
		* ```-op_gen``` the filename to save the consensus genotype matrix.

5. To evaluate the sensitivity, specificity use the following script:
	``` python ../evaluateMetrics.py -cg consensus_genotype.csv -gtG groundTruthFile -header true > eval_metrics.txt ```
For doublets run the script with a ```doublet``` flag:
	``` python ../evaluateMetrics.py -cg consensus_genotype.csv -gtG groundTruthFile -header true -doublet true -doubletFile doubletFileInformation > eval_metrics.txt ```

6. To get the V-measure use the following script:
	``` python ../evaluation.py -i "scite:"assignment.txt -G groundTruthFile -v >> eval_metrics.txt ```
For doublets run the script with a ```doublet``` flag:
	``` python ../evaluation.py -i "scite:"assignment.txt -G groundTruthFile -d true -df doubletFileInformation -v >> eval_metrics.txt ```

Look for example doublet file in the 'SCITE/example_data'

## <a name="sbmclone"></a>SBMClone ##

### Installing the software and preparing the files. These steps are general for all datasets, but they should be done only once. ###

1. Download SBMClone from https://github.com/raphael-group/SBMClone.git and the authors suggest using a conda environment to use it:
	``` 
	conda create -n sbmclone
	conda activate sbmclone
	conda install numpy scipy matplotlib networkx
	conda install -c conda-forge graph-tool 
	```
2. The input to SBMClone is a binary mutation matrix where rows correspond to cells, columns correspond to mutations, and each entry is a 1 if the corresponding cell has the corresponding mutation or 0 otherwise (equivalently, 0-entries could be represented as ?). The input format is a comma-separated text file in which each line encodes the row and column indices of a single 1-entry. Take a look at the example data in 'sbmclone/example_data/example.sbm.csv' to understand the input. 

3. Once SBMClone is installed you can run the script in the following way:
	* ``` python sbmclone.py example.sbm.csv -o outputDir ```
	* Parameters:
		* Pass the input file after the python script
		* ``` -o ```Pass the output directory to save the results 

4. SBMClone outputs the cells assigned to each clusters and we modify it such that we can use it for our evaluation script. To do that run the following script:
	* ``` python processOutput.py -input cluster-assignments.txt -output assignment.txt ```
	* Parameters:
		* ```-input``` Pass the output file by SBMClone.
		* ```-output``` Pass the output file name to save the results.

5. To get the V-measure use the following script:
	``` python ../evaluation.py -i "scclone:"assignment.txt -G groundTruthFile -v ```
We use "scclone" here because the assignment.txt is same as SCClone's assignment.txt file.

## <a name="ucsim"></a>UltraLow coverage simulator ##

### Steps to run the ultralow coverage simulator using different variables. ###

1. We use an existing tree as an input to the simulator to get the ultralow coverage datasets. There is an example tree in 'ultraLowCoverage_simulator/example_tree.csv' that can be used as an input.
2. Run the simulator using the following script:
	* ``` python sim_par.py -f example_tree.csv -c 4000 -n 5000 -cov 0.01 -P output_dir ```
	* Parameters:
		* ```-f``` pass the tree as input.
		* ```-c``` number of cells.
		* ```-n``` number of mutations.
		* ```-cov``` the coverage desired for the data.
		* ```-P``` the prefix and directory for output files.

## <a name="simulator"></a>Simulator ##

### Steps to run the simulator using different variables. ###

1. To generate the tree with beta-splitting variable and desired number of leaves use the following script:
	* ``` python gen_tree.py -F $i -B 0.2 -o output_tree.csv ```
	* Parameters:
		* ```-F``` pass the number of leaves or desired number of clones.
		* ```-B``` pass the beta split variable value.
		* ```-o``` output tree name to save.

2. Once we have the tree we can assign cells, mutations and include false positives, false negatives, missing data, doublets using the following script:
	* ``` python sim_par.py -a 0.01 -b 0.2 -m 0.2 -c 500 -n 200 -e 0.1 -f output_tree.csv -P outputFilePrefix ```
	* Parameters:
		* ```-a``` false positive rate.
		* ```-b``` false negative rate.
		* ```-m``` missing rate.
		* ```-c``` number of cells.
		* ```-n``` number of mutations.
		* ```-e``` doublet rate.
		* ```-f``` the tree generated earlier using number of leaves and beta splitting variable.
		* ```-P``` the prefix and directory for output files.



