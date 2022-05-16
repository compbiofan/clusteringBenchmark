import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input D matrix")
args = parser.parse_args()

input_df = pd.read_csv(args.input, sep='\t', index_col=0)
transposed_df = input_df.T

print(transposed_df)

# Save input for BnpC 
op_fname = (args.input.split('.'))[0]
transposed_df.to_csv(op_fname+'_bnpc.tsv', sep='\t', index=False, header=False)
