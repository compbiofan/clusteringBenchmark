import pandas as pd
import argparse

def convertToPandas(fileName,sim):
    pd_fromFile = pd.read_csv(fileName, sep='\t')
    if sim == "false":
        pd_fromFile = pd_fromFile.iloc[: , 1:] # Dropping the first column
    return pd_fromFile

def getClusterToCell(ca_file):
    clusterToCellDict = {}
    # First create a cluster to cell dictionary.
    with open(ca_file) as inp: # exp2_2_new.cell_assignment goes here.
        cell_clones_list = list(zip(*(line.strip().split(' ') for line in inp)))
    print(len(cell_clones_list))

    cell_no = 0
    for clones in cell_clones_list:
        clusterNo = clones[0]
        if clusterNo in clusterToCellDict:
            temp_list = clusterToCellDict[clusterNo]
            temp_list.append(cell_no)
            clusterToCellDict[clusterNo] = temp_list
        else:
            cell_list = []
            cell_list.append(cell_no)
            clusterToCellDict[clusterNo] = cell_list
        cell_no = cell_no + 1
    print(clusterToCellDict)
    return clusterToCellDict

def getGMatrix(cg_file,ca_file,opFile):
    clone_genotype_dict = {}
    with open(cg_file) as fp:
        Lines = fp.readlines()
        count = 0
        for line in Lines:
            genotype = (line.strip()).split('\t')
            clone_genotype_dict[count] = genotype
            count = count+1
    print("CG Dict ",clone_genotype_dict)

    with open(ca_file) as inp: # exp2_2_new.cell_assignment goes here.
        cell_clones_list = list(zip(*(line.strip().split(' ') for line in inp)))

    print(" Cell clones list ",cell_clones_list)
    GDoublePrimeDf = pd.DataFrame()
    count = 0
    for clone in cell_clones_list:
        GDoublePrimeDf[count] = clone_genotype_dict[int(clone[0])]
        count = count+1

    GDoublePrimeDf_T = GDoublePrimeDf.T
    print(GDoublePrimeDf_T)
    GDoublePrimeDf_T.to_csv(opFile, sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument("-cg", "--cg",dest ="cg", help="Input scclone_clone_genotypes file")
parser.add_argument("-ca", "--ca",dest ="ca", help="Input scclone_cell_assignment file")
parser.add_argument("-op", "-op",dest ="op", help="Ouput file to save")
args = parser.parse_args()

getGMatrix(args.cg,args.ca,args.op)

