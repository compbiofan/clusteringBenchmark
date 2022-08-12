import pandas as pd
import argparse
import numpy as np

def get_doublet_cells(doublet_info_file):
    #doublet_info_file = "sim_input/"+simDataType+"/"+repNo+"/input_"+simDataType+"_"+repNo+".SNVcell.csv"
    with open(doublet_info_file,'r') as fp:
        lines = fp.readlines()
        #print(lines)

    doublet_cells_list = []
    for line in lines:
        first_col = line.split('\t')
        if ';' in first_col[0]:
            cells = first_col[1]
            if ';' in cells:
                cells_list = (cells.strip()).split(';')
                doublet_cells_list.extend(cells_list)
            else:
                doublet_cells_list.append(cells.strip())
    print(" Doublet cells ",doublet_cells_list)
    return doublet_cells_list

def doublet_GT_CG(gt_df, cg_df, doublet_cells_list):
    cols = [i for i in range(0,len(gt_df.columns))]
    gt_df.columns = cols
    cg_df.columns = cols
    cg_df.index = gt_df.index
    print(" GT matrix ")
    print(gt_df)
    print(" CG matrix ")
    print(cg_df)
    #mutation_df = cluster_mutation_df.reindex(pos, columns=cells)
    for dCells in doublet_cells_list:
        gt_df.drop(int(dCells), inplace=True, axis=0)
        cg_df.drop(int(dCells), inplace=True, axis=0)
    #print(gt_df)
    #return gt_df

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
    print(CG_arr)
    print(" ============ ")
    print(GT_arr)
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

    print(" All correct entries ",correctNos," total 1s in GT ",total_1s," total 0s in GT ",total_0s)
    print(" All correct 1s ",correct1s)
    print(" All correct 0s ",correct0s)
    print(" All entries ",total_entries)
    accuracy = (correctNos/total_entries)*100
    sensitivity = (correct1s/total_1s)*100
    specificity = (correct0s/total_0s)*100
    print(" Accuracy ",accuracy)
    print(" Sensitivity ",sensitivity)
    print(" Specificity ",specificity)

parser = argparse.ArgumentParser()
parser.add_argument("-cg", "--cg",dest ="cg", help="Consensus genotype matrix")
parser.add_argument("-gtG", "--gtG",dest ="gtG", help="Ground truth matrix")
parser.add_argument("-sim", "--sim",dest ="sim", help="Simulated data")
parser.add_argument("-doublet", "--doublet",dest ="doublet", help="Doublet data")
parser.add_argument("-doubletFile", "--doubletFile",dest ="doubletFile", help="Doublet info file")
args = parser.parse_args()

if args.sim == "true":
    CG_df = sim_read_file(args.cg,"false")
    GT_df = sim_read_file(args.gtG,"true")
else:
    CG_df = read_file(args.cg)
    GT_df = read_file(args.gtG)

if args.doublet == "true":
    doublet_cells_list = get_doublet_cells(args.doubletFile)
    doublet_GT_CG(GT_df, CG_df, doublet_cells_list)

print(" CG =============================================== ")
print(CG_df)
print(" GT =============================================== ")
print(GT_df)
#calculate_metrics(CG_df, GT_df)
print(" ==================================================== ")
calculate_metric_from_vectors(CG_df,GT_df)
