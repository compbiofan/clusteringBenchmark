import pandas as pd
import argparse
import gzip
import shutil

# Get the cols based on the no. of mutations 
def getColumns():
    col_list = []
    for i in range(0,200):
        col_list.append('c'+str(i))
    return col_list

def process_input(input_file):
    input_df = pd.read_csv(args.input, sep='\t', header=None)
    input_df.index.name = 'cell_id'
    index_list = input_df.index.tolist()
    col_list = getColumns()
    input_df.columns = col_list
    print(input_df)
    return input_df

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input D matrix")
parser.add_argument("-output", "--output",dest ="output", help="Output to save the results")
args = parser.parse_args()

processed_df = process_input(args.input)

# Save input for BnpC
op_fname = (((args.input.split('/'))[4]).split('.'))[0]
print(op_fname)
processed_df.to_csv(args.output+'/'+op_fname+'.D.tsv', sep='\t')
with open(args.output+'/'+op_fname+'.D.tsv','rb') as f_in:
    with gzip.open(args.output+'/'+op_fname+'.D.tsv.gz','wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
