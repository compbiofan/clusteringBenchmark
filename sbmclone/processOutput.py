import argparse

def process_input(input_file, output_file):
    file=open(input_file,"r")
    line=file.readline().rstrip("\n")
    line=line.replace(","," ")
    op_file=open(output_file,"w")
    op_file.write(line)
    op_file.close()

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="cluster assignment file")
parser.add_argument("-output", "--output",dest ="output", help="Output to save the assignment file for V-measure")
args = parser.parse_args()
process_input(args.input, args.output)
