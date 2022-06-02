import pandas as pd
import argparse

def convertToPandas(fileName,sim):
    pd_fromFile = pd.read_csv(fileName, sep='\t')
    if sim == "false":
        pd_fromFile = pd_fromFile.iloc[: , 1:] # Dropping the first column
    return pd_fromFile

def getGMatrix(inputFile,cg_file,ca_file,opFile):
    original_df = pd.read_csv(inputFile, sep='\t', index_col = 0)
    clone_genotype_df = convertToPandas(cg_file,"false") # exp2_2_new.clone_genotypes goes here.
    print(" Clone genotype ")
    print(clone_genotype_df)
    with open(ca_file) as inp: # exp2_2_new.cell_assignment goes here.
        cell_clones_list = list(zip(*(line.strip().split(' ') for line in inp)))

    col_list = list(original_df.columns)
    cell_clones_list = cell_clones_list[1:]
    print(" Cell clones list ",len(cell_clones_list))
    print("Index ",original_df.index)
    print("Columns ",col_list)

    GDoublePrimeDf = pd.DataFrame()
    count = 0
    for clone in cell_clones_list:
        df_index = int(clone[0]) - 1
        #print(clone," ",clone_genotype_df.iloc[df_index])
        GDoublePrimeDf[count] = clone_genotype_df.iloc[df_index]
        count = count+1

    GDoublePrimeDf_T = GDoublePrimeDf.T
    GDoublePrimeDf_T.columns = list(original_df.columns)
    GDoublePrimeDf_T.index = list(original_df.index)
    #print(GDoublePrimeDf_T)
    GDoublePrimeDf_T.to_csv(opFile, sep='\t') # output file to save the G matrix

def sim_getGMatrix(cg_file,ca_file,opFile):
    clone_genotype_df = convertToPandas(cg_file,"true") # exp2_2_new.clone_genotypes goes here.
    print(" Clone genotype ")
    with open(ca_file) as inp: # exp2_2_new.cell_assignment goes here.
        cell_clones_list = list(zip(*(line.strip().split(' ') for line in inp)))

    print(" Cell clones list ",len(cell_clones_list))
    GDoublePrimeDf = pd.DataFrame()
    count = 0
    for clone in cell_clones_list:
        #print(" Count ",count)
        #if count >= len(col_list):
        #    continue
        #print(col_list[count])
        df_index = int(clone[0]) - 1
        #print(clone," ",clone_genotype_df.iloc[df_index])
        GDoublePrimeDf[count] = clone_genotype_df.iloc[df_index]
        count = count+1

    GDoublePrimeDf_T = GDoublePrimeDf.T
    #GDoublePrimeDf_T.index = list(original_df.index)
    print(GDoublePrimeDf_T)
    GDoublePrimeDf_T.to_csv(opFile, sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input D file")
parser.add_argument("-cg", "--cg",dest ="cg", help="Input scclone_clone_genotypes file")
parser.add_argument("-ca", "--ca",dest ="ca", help="Input scclone_cell_assignment file")
parser.add_argument("-sim", "--sim",dest ="sim", help="Simulated data")
parser.add_argument("-op", "-op",dest ="op", help="Ouput file to save")
args = parser.parse_args()

if args.sim == "true":
    sim_getGMatrix(args.cg,args.ca,args.op)
else:
    getGMatrix(args.input,args.cg,args.ca,args.op)
