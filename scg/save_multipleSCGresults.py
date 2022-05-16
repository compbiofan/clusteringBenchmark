import yaml
import subprocess
import argparse
import pandas as pd
import os
import pathlib
from pathlib import Path
import random
import shutil
from shutil import rmtree

from glob import glob

def read_yaml(yaml_file):
    with open(yaml_file,'r') as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        config_params = yaml.load(file, Loader=yaml.FullLoader)
        print(config_params['num_clusters'])
        return config_params
        #print(config_params['num_clusters'])

def calculate_max_clusters(fileName):
#    cluster_list = []
    df_name = pd.read_csv(fileName, sep='\t', index_col = 0)
    no_of_cells = len(list(df_name.index.values))
#    print(" No of cells ",no_of_cells)
    max_clusters = int(no_of_cells/4)
#    start_cluster = int(max_clusters/4)
#    print(" Initial cluster no ",start_cluster)
#    cluster_list.append(start_cluster)
#    initial_sc = start_cluster
#    for i in range(start_cluster,max_clusters,initial_sc):
#        print(" Start cluster ",i)
#        start_cluster = int(start_cluster*2)
#        #print("Start cluster ",start_cluster)
#        cluster_list.append(start_cluster)

#    print(" Cluster list ",cluster_list)
    return max_clusters

# Save the lower bound values and seed value in a dictionary
def save_seed_lb(opDir, lb_dir, seed_value, seed_lb_dict):
    path = opDir+'/'+lb_dir
    print(path)
    #os.chdir(path)
    for file in os.listdir(path):
        if file.endswith(".txt"):
            file_path = path+"/"+file
            with open(file_path, 'r') as f:
                lines = [line.rstrip() for line in f]
                converge = (lines[0].split(':'))[1].strip()
                lower_bound = (lines[1].split(':'))[1].strip()
                if converge == "true":
                    seed_lb_dict[seed_value] = lower_bound
    return seed_lb_dict

def get_key(elbo,seed_lb_dict):
    for key, value in seed_lb_dict.items():
         if elbo == float(value):
             return key

def find_best_run(config_params, fileName, cluster, nIters, opDir, seed_lb_dict,config_file_name):
    best_seed = 0
    elbo = float(list(seed_lb_dict.values())[0])
    for key,value in seed_lb_dict.items():
        if float(value) > elbo:
            elbo = float(value)
            best_seed = key
        if best_seed == 0:
            best_seed = get_key(elbo,seed_lb_dict)
    print("Max elbo ",elbo," for seed ",best_seed)
    # After this remove existing directories and run scg again with new seed value.
    pattern = os.path.join(opDir, "scg_clusterNo_*")
    for item in glob(pattern):
        if not os.path.isdir(item):
            continue
        rmtree(item)
    scg_cmd = 'time scg run_singlet_model --config_file ../../SCG_Roth/scg/examples/config_files/'+config_file_name+' --seed '+str(best_seed)+' --lower_bound_file '+opDir+'/'+'lower_bound.txt --max_num_iters '+nIters+' --out_dir '+opDir
    subprocess.call(scg_cmd, shell=True)
        
def save_multipleSCGresults(config_params, fileName, cluster, nIters, opDir):
    #for i in clusterNo_list:
    cluster_update = {'num_clusters': int(cluster)}
    #ufileName = {'data': {'snv': {'file': fileName, 'gamma_prior':[[98, 1, 1], [25, 50, 25], [1, 1, 98]], 'state_prior': [1, 1, 1]}}}
    ufileName = {'data': {'snv': {'file': fileName, 'gamma_prior':[[9.99, 0.01, 1.0e-15], [2.5, 7.5, 1.0e-15], [1.0e-15, 1.0e-15, 1]], 'state_prior': [1, 1, 1]}}}
    #gamma_prior = {'data': {'snv': {'gamma_prior':[[98, 1, 1], [25, 50, 25], [1, 1, 98]]}}}
    #state_prior = {'data': {'snv': {'state_prior': [1, 1, 1]}}}
    config_name = (((fileName.split('/'))[2]).split('.'))[0]
    config_file_name = 'config_'+config_name+'.yaml'
    print(" Config name ",config_file_name)
    config_params.update(cluster_update)
    config_params.update(ufileName)
    #print(config_params)
    # Update the config file to have the values
    with open('../../SCG_Roth/scg/examples/config_files/'+config_file_name, 'w') as fp:
        documents = yaml.dump(config_params, fp)
    
    seed_lb_dict = {}
    for j in range(0, 4): # Based on the restart value change this loop. Otherwise let it stay as it is.
        print("Iteration ",j)
        seed_value = random.randint(0, 10000)
        print(seed_value)
        dir_name = 'scg_clusterNo_'+str(cluster)+'_iter_'+str(j)
        #dir_name = 'scg_clusterNo_'+str(cluster)
        #wd = os.getcwd()
        #print(" Current working dir ",wd," change to ",opDir)
        #os.chdir(opDir)
        print(" Dir name ",opDir+'/'+dir_name)
        subprocess.call('mkdir '+opDir+'/'+dir_name, shell=True)
        scg_cmd = 'time scg run_singlet_model --config_file ../../SCG_Roth/scg/examples/config_files/'+config_file_name+' --seed '+str(seed_value)+' --lower_bound_file '+opDir+'/'+dir_name+'/lower_bound_file'+str(seed_value)+'.txt --max_num_iters '+nIters+' --out_dir '+opDir+'/'+dir_name
        subprocess.call(scg_cmd, shell=True)
        #os.chdir(wd)
        seed_lb_dict = save_seed_lb(opDir, dir_name, seed_value, seed_lb_dict)
    print(seed_lb_dict)
    find_best_run(config_params, fileName, cluster, nIters, opDir, seed_lb_dict, config_file_name)

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (similar to input to SCG)")
#parser.add_argument("-cluster", "--cluster",nargs='+', dest = "cluster_list", help="List of cluster numbers")
parser.add_argument("-scg_config","--scg_config",dest="scg_config", help="SCG Config yaml file")
parser.add_argument("-niters", "--niters",dest ="niters", help="No of iterations")
parser.add_argument("-opDir", "--opDir",dest ="opDir", help="Output Directory")
#args = parser.parse_args()
args, unknown = parser.parse_known_args()

cluster = calculate_max_clusters(args.input)

config_params = read_yaml(args.scg_config)
#clusterNo_list = args.cluster_list
save_multipleSCGresults(config_params, args.input, cluster, args.niters, args.opDir)
 
#config_params = read_yaml('../../SCG_Roth/scg/examples/config.yaml')
#clusterNo_list = [5,6,7,8,9,10]
#save_multipleSCGresults(config_params, '/gpfs/research/fangroup/rk18g/longitudinal/LACE-UTILITIES/real_data/Breast_Cancer/p494/final_LACE_input_data/merged.tsv.gz',clusterNo_list, 10)

# update dictionary config_params and dump to the same yaml file
