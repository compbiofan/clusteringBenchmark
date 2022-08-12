import pandas as pd
import argparse

def convertToPandas(fileName,sim):
    pd_fromFile = pd.read_csv(fileName, sep='\t')
    if sim == "false":
        pd_fromFile = pd_fromFile.iloc[: , 1:] # Dropping the first column
    return pd_fromFile

#''' Vote to break the tie if there are same probabilities for the genotypes. '''
#def voting_gt(cell_gt_list,pos):
#    count1 = 0
#    count0 = 0
#    count3 = 0 # 3 is for missing data
    #print(cell_gt_list)
    # pos is the mutation index till 200 and an outer loop before the below for loop
#    for cell in cell_gt_list:
        #print(cell," ",pos," ",cell[0][pos])
#        if cell[0][pos] == 1:
            #print(" Condn 1 satisfied ")
#            count1 = count1+1
#        elif cell[0][pos] == 0:
#            count0 = count0+1
#        elif cell[0][pos] == 3:
            #print(" Condn 3 satisfied ")
#            count3 = count3+1
    #print(" Count1 ",count1," Count 0 ",count0," Count 3 ",count3)
#    if count1 > count0 or count1 == count0 or count1 > count3: # What if there is a tie with same no of 1s and 0s? Select randomly.
#        return 1
#    elif count0 > count1 or count0 > count3:
#        return 0
#    elif count3 > count1 or count3 > count0:
#        return 0

#''' Return the genotype based on voting. '''
#def breakTie_WithMajority(clusterToCellDict,D_table):
    # From the cell To cluster dict get the cell IDs for each cluster. Then look for the genotypes of those cell IDs.
    # Based on the pos select the majority 
    #print(clusterToCellDict)
#    cluster_genotype = {}
#    for cluster, cells in clusterToCellDict.items():
#        cell_genotype_list = []
#        for cell in cells:
            #print(cell)
            # Get the rows of D based on indexID cell and append to a list. From that list for the given pos find the majority of the column.
#            g_list = list(D_table.iloc[[cell]].values)
            #print(g_list)
#            cell_genotype_list.append(g_list)
        #print(cluster," ",cell_genotype_list)
#        voted_gt_list = []
#        for pos in range(0,200):
#            voted_gt = voting_gt(cell_genotype_list,pos)
#            voted_gt_list.append(voted_gt)
        #print(cluster," ",voted_gt_list)
#        cluster_genotype[cluster] = voted_gt_list
#    print(cluster_genotype)
    #return voted_gt

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
    #clone_genotype_df = convertToPandas(cg_file,"true") # exp2_2_new.clone_genotypes goes here.
    #print(" Clone genotype ")
    
    #print(" Clone genotype DF ======== ")
    #print(clone_genotype_df)

    clone_genotype_dict = {}
    with open(cg_file) as fp:
        Lines = fp.readlines()
        count = 0
        for line in Lines:
            genotype = (line.strip()).split('\t')
            clone_genotype_dict[count] = genotype
            count = count+1
    print("CG Dict ",clone_genotype_dict)

    #cg_df_index = list(clone_genotype_df.index)
    #print(list(clone_genotype_df.index))
    with open(ca_file) as inp: # exp2_2_new.cell_assignment goes here.
        cell_clones_list = list(zip(*(line.strip().split(' ') for line in inp)))

    print(" Cell clones list ",cell_clones_list)
    GDoublePrimeDf = pd.DataFrame()
    count = 0
    for clone in cell_clones_list:
        #print(" Clone ",clone)
        #print(" Count ",count)
        #if count >= len(col_list):
        #    continue
        #print(col_list[count])

        #df_index = int(clone[0]) - 1
        
        #if df_index not in cg_df_index:
        #    continue
        #print(" DF index ",df_index)
        #print(clone," ",clone_genotype_df.iloc[df_index])

        #GDoublePrimeDf[count] = clone_genotype_df.iloc[df_index]
        GDoublePrimeDf[count] = clone_genotype_dict[int(clone[0])]
        count = count+1

    #print(" GDoublePrimeDf matrix ====== ")
    #print(GDoublePrimeDf)
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

#D_table = pd.read_csv(args.input, sep='\t', header=None)
#print(D_table)

#clusterToCell_dict = getClusterToCell(args.ca)
#breakTie_WithMajority(clusterToCell_dict,D_table)

if args.sim == "true":
    sim_getGMatrix(args.cg,args.ca,args.op)
else:
    getGMatrix(args.input,args.cg,args.ca,args.op)
