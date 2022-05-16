import pandas as pd
import argparse

def convertToPandas(fileName):
    pd_fromFile = pd.read_csv(fileName, sep='\t')
    pd_fromFile = pd_fromFile.iloc[: , 1:] # Dropping the first column
    return pd_fromFile

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input D file")
parser.add_argument("-cg", "--cg",dest ="cg", help="Input scclone_clone_genotypes file")
parser.add_argument("-ca", "--ca",dest ="ca", help="Input scclone_cell_assignment file")
parser.add_argument("-op", "-op",dest ="op", help="Ouput file to save")
args = parser.parse_args()

original_df = pd.read_csv(args.input, sep='\t', index_col = 0)
clone_genotype_df = convertToPandas(args.cg) # exp2_2_new.clone_genotypes goes here.
print(" Clone genotype ")
print(clone_genotype_df)

#fname = open("exp2_1.cell_assignment")
#for l in fname.readlines():
#    cell_clones_str = l.strip().split("\t")
#fname.close()

#cell_clones_list = cell_clones_str[0].split("\t")
#print(cell_clones_list)

with open(args.ca) as inp: # exp2_2_new.cell_assignment goes here.
    cell_clones_list = list(zip(*(line.strip().split(' ') for line in inp))) 

col_list = list(original_df.columns)
cell_clones_list = cell_clones_list[1:]
print(" Cell clones list ",len(cell_clones_list))
print("Index ",original_df.index)
print("Columns ",col_list)

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

#print(GDoublePrimeDf.T.columns)
#GDoublePrimeDf_T = pd.DataFrame(columns = original_df.columns, index = original_df.index)
GDoublePrimeDf_T = GDoublePrimeDf.T
GDoublePrimeDf_T.columns = list(original_df.columns)
GDoublePrimeDf_T.index = list(original_df.index)
#print(GDoublePrimeDf_T)
GDoublePrimeDf_T.to_csv(args.op, sep='\t') # output file to save the G matrix

#GDoublePrimeDf_T_1 = GDoublePrimeDf_T.reindex(original_df.index, columns = original_df.columns)
#print(GDoublePrimeDf_T_1)
