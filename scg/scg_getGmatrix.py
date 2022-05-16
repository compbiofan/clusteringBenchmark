import pandas as pd
import json
from collections import defaultdict
import argparse

def convertToPandas(fileName):
    pd_fromFile = pd.read_csv(fileName, sep='\t', index_col=0)
    return pd_fromFile

def cell_cluster(cluster_posterior):
    tsv_data = pd.read_csv(cluster_posterior, compression='gzip', sep='\t', index_col=0)
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

''' Function to get the G' matrix. Get it from SCG genotype_posterior and replace SNV 2 with 1.'''
def build_GprimeMatrix(gp_table, cellToClusterDict):
    clusterToCellDict = {} # Building a dictionary with cluster IDs as keys and cell IDs as correspondig values
    #print(" GP table ",gp_table)
    #print(" cellToClusterDict ",cellToClusterDict)
    for cellID in cellToClusterDict:
        clusterID = cellToClusterDict[cellID]
        if clusterID in clusterToCellDict:
            temp_list = clusterToCellDict[clusterID];
            temp_list.append(cellID)
        else:
            clusterToCellDict[clusterID] = [cellID]
    #print("Cluster to cell Dict ",clusterToCellDict)

    GprimeDf = pd.DataFrame()
    for cluster in clusterToCellDict:
        for row in gp_table.itertuples():
            if cluster == str(row.Index) and getattr(row,'event_type') == 'snv' and getattr(row,'probability') > 0.80:
                #prob = getattr(row,'probability')
                #if math.ceil(prob) == 1:
                    pos = getattr(row,'event_id')
                    SNV = getattr(row,'event_value')
                    if SNV == 2.0:
                        SNV = 1.0
                    GprimeDf.loc[cluster,pos] = SNV

    return GprimeDf

''' Function to build the G'' matrix from the output of SCG. '''
def build_GDoublePrimeMatrix(gp_table, cellToClusterDict):
    gDoublePrimeMatrix = pd.DataFrame()
    for cell in cellToClusterDict: # key is the Cell ID so we check for each cell ID its SNV types on a specific position
        clusterID = cellToClusterDict[cell]
        for row in gp_table.itertuples():
            if clusterID == str(row.Index) and getattr(row,'event_type') == 'snv' and getattr(row,'probability') > 0.80:
                #prob = getattr(row,'probability')
                #if math.ceil(prob) == 1:
                    pos = getattr(row,'event_id')
                    SNV = getattr(row,'event_value')
                    gDoublePrimeMatrix.loc[cell,pos] = int(SNV)
        #for index, row in gp_table.iterrows():
        #    if clusterID == str(index) and row['event_type'] == 'snv': # Choose the cluster the cell ID belongs to and get the position(event_id) and SNV(event_value) to fill the matrix
        #        pos = row['event_id']
        #        SNV = row['event_value']
                #print("Event value .. ",event_value)
        #        initial_matrix.loc[cell,pos] = SNV
    return gDoublePrimeMatrix

parser = argparse.ArgumentParser()
parser.add_argument("-cp", "--cluster",dest ="cluster", help="Cluster posterior filename.")
parser.add_argument("-gp", "--genotype",dest ="genotype", help="Genotype posterior from SCG")
parser.add_argument("-output", "--output",dest ="output", help="Consensus genotype output filename")
args = parser.parse_args()
cellToClusterDict = cell_cluster(args.cluster)

gp_path = args.genotype
gp_table = convertToPandas(gp_path) # gp_table is the dataframe with values from genotype_posteriors.tsv

cp_path = args.cluster
cp_table = convertToPandas(cp_path) # cp_table is the dataframe with values from cluster_posteriors.tsv
cp_table['cell_id'] = cp_table.index # cell IDs are the index here

GDoublePrimeDf = build_GDoublePrimeMatrix(gp_table, cellToClusterDict)
print(GDoublePrimeDf)
GDoublePrimeDf.to_csv(args.output, sep='\t')
