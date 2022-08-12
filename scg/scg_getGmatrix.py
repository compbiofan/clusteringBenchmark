import pandas as pd
import json
from collections import defaultdict
import argparse

def convertToPandas(fileName):
    pd_fromFile = pd.read_csv(fileName, sep='\t', index_col=0)
    return pd_fromFile

''' Get the final consensus genotype matrix using the method parameters.'''
def get_CGMatrix(cellToClusterDict, D_table, cluster_genotype):
    #print(D_table.columns)
    #print("No of cells ",len(cellToClusterDict))
    gDoublePrimeMatrix = pd.DataFrame(index=cellToClusterDict.keys(),columns=D_table.columns)
    for cell in cellToClusterDict:
        clusterID = cellToClusterDict[cell]
        genotype = cluster_genotype[clusterID]
        gDoublePrimeMatrix.iloc[cell] = genotype
    return gDoublePrimeMatrix

''' Get the dictionary where cluster ID is the key and cells belonging to that cluster is the value. '''
def get_ClusterToCell(cellToClusterDict):
    clusterToCellDict = {}
    for key,value in cellToClusterDict.items():
        if value in clusterToCellDict:
            temp_list = clusterToCellDict[value]
            temp_list.append(key)
            clusterToCellDict[value] = temp_list
        else:
            cell_list = []
            cell_list.append(key)
            clusterToCellDict[value] = cell_list
    print(clusterToCellDict.keys())
    return clusterToCellDict

''' Vote to break the tie if there are same probabilities for the genotypes. '''
def voting_gt(cell_gt_list,pos):
    count1 = 0
    count0 = 0
    count3 = 0 # 3 is for missing data
    #print(cell_gt_list)
    for cell in cell_gt_list:
        #print(cell," ",pos," ",cell[0][pos])
        if cell[0][pos] == 1:
            #print(" Condn 1 satisfied ")
            count1 = count1+1
        elif cell[0][pos] == 0:
            count0 = count0+1
        elif cell[0][pos] == 3:
            #print(" Condn 3 satisfied ")
            count3 = count3+1
    #print(" Count1 ",count1," Count 0 ",count0," Count 3 ",count3)
    if count1 > count0 or count1 > count3: # What if there is a tie with same no of 1s and 0s? Select randomly.
        return 1
    elif count0 > count1 or count0 > count3:
        return 0
    elif count1 == count0 or count3 > count1 or count3 > count0:
        return "Random"

''' Return the genotype based on voting. '''
def breakTie_WithMajority(clusterToCellDict,clusterId,D_table,pos):
    # From the cell To cluster dict get the cell IDs for each cluster. Then look for the genotypes of those cell IDs.
    # Based on the pos select the majority 
    #print(clusterToCellDict)
    cells = clusterToCellDict[clusterId]
    cell_genotype_list = []
    for cell in cells:
        # Get the rows of D based on indexID cell and append to a list. From that list for the given pos find the majority of the column.
        g_list = list(D_table.iloc[[cell]].values)
        cell_genotype_list.append(g_list)
    #print(clusterId," ",cell_genotype_list)
    voted_gt = voting_gt(cell_genotype_list,pos)
    return voted_gt 

''' From the genotype posterior file select the correct genotypes for the clusters. '''
def getG_fromGenotypePosterior(gp_table,cellToClusterDict,D_table):
    # if event_type == 'snv' then create a dictionary with cluster ID as the key and value as lists of 'pos-genotype-prob'.
    # Then access this dictionary to get the lists and hover over to get the genotypes of the probabilities based on highest probability values.
    gp_dict = {}
    for row in gp_table.itertuples():
        if getattr(row,'event_type') == 'snv':
            key = str(row.Index)+getattr(row,'event_id')
            #print(type(getattr(row,'probability')))
            if key in gp_dict:
                temp_list = gp_dict.get(key)
                temp_list.append(str(getattr(row,'event_value'))+':'+str(getattr(row,'probability')))
                gp_dict[key] = temp_list
            else:
                list_val = []
                list_val.append(str(getattr(row,'event_value'))+':'+str(getattr(row,'probability')))
                gp_dict[key] = list_val
    #print(gp_dict.keys())
    cluster_genotype = {}
    #genotype_list = [0 for x in range(200)]
    clusterToCellDict = get_ClusterToCell(cellToClusterDict)
    for clusterPos,genotype_prob in gp_dict.items():
        clusterId = (clusterPos.split('c'))[0]
        if clusterId not in clusterToCellDict.keys():
            continue
        pos = int((clusterPos.split('c'))[1])

        prob_list = []
        for gprob in genotype_prob:
            genotype = (gprob.split(':'))[0]
            prob = (gprob.split(':'))[1]
            prob_list.append(float(prob))

        max_prob = max(prob_list)
        max_prob_ties = [i for i, j in enumerate(prob_list) if j == max_prob]
        #print(" Max prob ",max_prob," Max prob ties ",max_prob_ties)
        if len(max_prob_ties) > 1:
            #print(" Cluster ID ",clusterId," for position ",pos," max prob ",max(prob_list))
            tie_g = breakTie_WithMajority(clusterToCellDict,clusterId,D_table,pos)
            #print(" Tie G ",tie_g)
            if tie_g == 0 or tie_g == 1:
                max_g = tie_g
                #print(" Cluster ID ",clusterId," for position ",pos," max g ",max_g)
            else:
                max_g = prob_list.index(max_prob)
        elif len(max_prob_ties) == 1:
            max_g = prob_list.index(max_prob)
        if clusterId in cluster_genotype:
            temp_list = []
            temp_list = cluster_genotype[clusterId]
            #print(temp_list)
            temp_list[pos] = max_g
            cluster_genotype[clusterId] = temp_list
        else:
            genotype_list = [0 for x in range(200)]
            genotype_list[pos] = max_g
            cluster_genotype[clusterId] = genotype_list
    return cluster_genotype

''' Function to assign cluster nos to cells. ''' 
def cell_cluster(cluster_posterior):
    tsv_data = pd.read_csv(cluster_posterior, sep='\t', index_col=0)
    #print(tsv_data)
    # get the column of max values in every row and make new column with max values
    #print(tsv_data.idxmax(axis=1))
    tsv_data['clusters'] = tsv_data.idxmax(axis=1)
    tsv_data['cell_id'] = tsv_data.index # converting cell_id index to column in dataframe
    cellToClusterDict = {}
    #print(tsv_data)
    for idx,row in tsv_data.iterrows():
        cellToClusterDict[row['cell_id']] = row['clusters']

    print(" No of clusters in the file ",len(set(cellToClusterDict.values())))
    return cellToClusterDict

parser = argparse.ArgumentParser()
parser.add_argument("-cp", "--cluster",dest ="cluster", help="Cluster posterior filename.")
parser.add_argument("-gp", "--genotype",dest ="genotype", help="Genotype posterior from SCG")
parser.add_argument("-D", "--D",dest ="D", help="Input matrix")
parser.add_argument("-output", "--output",dest ="output", help="Consensus genotype output filename")
args = parser.parse_args()
cellToClusterDict = cell_cluster(args.cluster)

gp_path = args.genotype
gp_table = convertToPandas(gp_path) # gp_table is the dataframe with values from genotype_posteriors.tsv
#print(gp_table)

D_table = convertToPandas(args.D)
#print(D_table)
#cp_path = args.cluster
#cp_table = convertToPandas(cp_path) # cp_table is the dataframe with values from cluster_posteriors.tsv
#cp_table['cell_id'] = cp_table.index # cell IDs are the index here

cluster_genotype_dict = getG_fromGenotypePosterior(gp_table,cellToClusterDict,D_table)
gDoublePrimeMatrix = get_CGMatrix(cellToClusterDict, D_table, cluster_genotype_dict)
print(gDoublePrimeMatrix)
print(" Output fileName ",args.output)
gDoublePrimeMatrix.to_csv(args.output, sep='\t')
