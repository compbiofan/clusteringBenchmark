# This file gets the G matrix with the cells and mutations.
import pandas as pd
import argparse

def get_cell_cluster(assignment):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cell_cluster = {}
    print(len(assignment_list))
    cell_count = 1
    for i in assignment_list:
        if '\n' in i:
            i = i.replace("\n","")
        cell_cluster[cell_count] = i
        cell_count = cell_count+1
    return cell_cluster

def get_cell_mutation(mutations,input_D):
    cluster_mutation_df = pd.read_csv(mutations, sep='\t', index_col = 0)
    #cluster_mutation_df_T = cluster_mutation_df.T
    print(cluster_mutation_df)
    input_df = pd.read_csv(input_D, sep='\t', index_col = 0)
    #print(input_df)
    cells = list(input_df.index)
    pos = list(input_df.columns)
    cluster_mutation_df.index = pos
    cluster_mutation_df.columns = cells
    #mutation_df = cluster_mutation_df.reindex(pos, columns=cells)
    print(cluster_mutation_df)
    return cluster_mutation_df.T

parser = argparse.ArgumentParser()
parser.add_argument("-cc", "--cc",dest ="cc", help="Assignment.txt file indicating clusters of cells")
parser.add_argument("-gp", "--gp",dest ="gp", help="Mutations for each cluster")
parser.add_argument("-input", "--input",dest="input", help="Sample input to bnpc")
parser.add_argument("-op","--op",dest="op", help="Output file to save")
args = parser.parse_args()

cell_cluster = get_cell_cluster(args.cc)
#print(cell_cluster.values())
cell_mutation_df = get_cell_mutation(args.gp,args.input)
cell_mutation_df.to_csv(args.op,sep='\t')


