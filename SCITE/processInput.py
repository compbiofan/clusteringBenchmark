import pandas as pd
import argparse

def process_input(input_file, output_file):
    file = open(input_file,"r")
    line = file.readline().rstrip('\n')
    modified_lines = []
    while(line != ""):
        line = line.replace("\t"," ")
        modified_lines.append(line)
        line = file.readline().rstrip('\n')
    with open(output_file,"w+") as f:
        for l in modified_lines:
            f.write(l)
            f.write("\n")

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input data matrix")
parser.add_argument("-output", "--output",dest ="output", help="Output data matrix")
args = parser.parse_args()

input_df = pd.read_csv(args.input, sep='\t', header=None)
transposed_df = input_df.T
transposed_df.to_csv(args.output, sep=' ', index=False, header=False)

#process_input(args.input, args.output)
