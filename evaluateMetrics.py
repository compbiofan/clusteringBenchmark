import pandas as pd
import argparse

def read_file(tsvFile):
    df = pd.read_csv(tsvFile, sep='\t', index_col=0)
    return df

def calculate_metrics(CG_df, GT_df):
    pos = list(CG_df.columns.values)
    cells = list(CG_df.index.values)
    #print(pos)
    #print(cells)
    correctNos = 0
    correct1s = 0
    correct0s = 0
    total_entries = len(pos) * len(cells)
    total_1s = 0
    total_0s = 0
    for p1 in pos:
        for row in GT_df.itertuples():
            cell_value = getattr(row,p1)
            if cell_value == 1:
                total_1s = total_1s+1
            if cell_value == 0:
                total_0s = total_0s+1
    for p1 in pos:
        for row in CG_df.itertuples():
            cg = getattr(row,p1)
            #print(" Int CG ",int(cg))
            gt = GT_df.at[row.Index,p1]
            if cg == gt:
                correctNos = correctNos+1
            if cg == 1 and gt == 1:
                correct1s = correct1s+1
            if cg == 0 and gt == 0:
                correct0s = correct0s+1

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
args = parser.parse_args()
CG_df = read_file(args.cg)
GT_df = read_file(args.gtG)
print(CG_df)
print(GT_df)
calculate_metrics(CG_df, GT_df)
