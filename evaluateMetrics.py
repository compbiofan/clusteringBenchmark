import pandas as pd
import argparse
import numpy as np

def sim_read_file(tsvFile,gt):
    if gt == "true":
        df = pd.read_csv(tsvFile, sep='\t', header=None)
    else:
        df = pd.read_csv(tsvFile, sep='\t', index_col=0)
    return df

def read_file(tsvFile):
    df = pd.read_csv(tsvFile, sep='\t', index_col=0)
    return df

def calculate_metric_from_vectors(CG_df,GT_df):
    CG_arr = CG_df.to_numpy()
    GT_arr = GT_df.to_numpy()
    correctNos = 0
    correct1s = 0
    correct0s = 0
    total_1s = 0
    total_0s = 0
    num_cells, num_mutation = CG_arr.shape
    total_entries = num_cells * num_mutation
    for x, y in np.ndindex(CG_arr.shape):
        if CG_arr[x,y] == GT_arr[x,y]:
            correctNos = correctNos+1
        if CG_arr[x,y] == 1 and GT_arr[x,y] == 1:
            correct1s = correct1s+1
        if CG_arr[x,y] == 0 and GT_arr[x,y] == 0:
            correct0s = correct0s+1

    for x, y in np.ndindex(GT_arr.shape):
        if GT_arr[x,y] == 1:
            total_1s = total_1s+1
        if GT_arr[x,y] == 0:
            total_0s = total_0s+1
    accuracy = (correctNos/total_entries)*100
    sensitivity = (correct1s/total_1s)*100
    specificity = (correct0s/total_0s)*100
    print(" All correct entries ",correctNos," total 1s in GT ",total_1s," total 0s in GT ",total_0s)
    print(" All correct 1s ",correct1s)
    print(" All correct 0s ",correct0s)
    print(" All entries ",total_entries)
    print(" Accuracy ",accuracy)
    print(" Sensitivity ",sensitivity)
    print(" Specificity ",specificity)

parser = argparse.ArgumentParser()
parser.add_argument("-cg", "--cg",dest ="cg", help="Consensus genotype matrix")
parser.add_argument("-gtG", "--gtG",dest ="gtG", help="Ground truth matrix")
parser.add_argument("-sim", "--sim",dest ="sim", help="Simulated data")
args = parser.parse_args()

if args.sim == "true":
    CG_df = sim_read_file(args.cg,"false")
    GT_df = sim_read_file(args.gtG,"true")
else:
    CG_df = read_file(args.cg)
    GT_df = read_file(args.gtG)

print(CG_df)
print(GT_df)
#calculate_metrics(CG_df, GT_df)
print(" ==================================================== ")
calculate_metric_from_vectors(CG_df,GT_df)
