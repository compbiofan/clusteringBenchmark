import argparse
import pandas as pd

''' Map the clones genotypes or the cells to each clones. '''
def map_clones(input_file):
    file = open(input_file,"r")
    line = file.readline().rstrip('\n')
    clone_no = 0
    clone_dict = {}
    while(line != ""):
       clone_dict[clone_no] = line
       clone_no = clone_no+1
       line = file.readline().rstrip('\n') 
    print(" Clone dict ",clone_dict)
    return clone_dict

''' Create an assignment file to calculate V-measure. '''
def create_assignment_file(clones_cells, noOfCells, op_ass):
    cell_assignment = [-1] * int(noOfCells)
    #print("=============================")
    for clones, cells in clones_cells.items():
        #print("Clone ",clones," cells ",cells)
        cell_list = cells.split(',')
        for cell in cell_list:
            #cell_assignment.insert(int(cell), clones)
            cell_assignment[int(cell)-1] = clones
        #print(cell_assignment)
    print(" Cell assignment ",cell_assignment)
    f = open(op_ass, "w")
    f.write(' '.join(list(map(str,cell_assignment))))
    f.close()

''' Get cells genotype from the clonal genotype. '''
def get_cells_genotype(clone_cells,clone_genotype,op_file):
    cells_genotype = {}
    for clone, cells in clone_cells.items():
        cells_list = cells.split(",")
        clonegen = (clone_genotype[clone]).split(" ")
        print(len(clone_genotype[clone]))
        for cell in cells_list:
            cells_genotype[int(cell)] = clonegen
    
    cell_keys = list(cells_genotype.keys())
    cell_keys.sort()
    sorted_dict = {i: cells_genotype[i] for i in cell_keys}
    print(sorted_dict.keys())
    output_df = pd.DataFrame.from_dict(sorted_dict, orient='index')
    print(output_df)
    output_df.to_csv(op_file,sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument("-cells", "--cells",dest ="cells", help="Cells assigned to each clones")
parser.add_argument("-cloneg", "--cloneg",dest ="cloneg", help="Clonal genotype")
parser.add_argument("-noOfCells", "--noOfCells",dest="noOfCells", help="No of cells")
parser.add_argument("-output", "--output",dest="output", help="Output file")
parser.add_argument("-op_ass", "--op_ass",dest="op_ass", help="Assignment file")
args = parser.parse_args()

clone_cells = map_clones(args.cells)
clone_genotype = map_clones(args.cloneg)
#print(" Cells of clones ",clone_cells)
#print(" Clone genotypes ",clone_genotype)
get_cells_genotype(clone_cells,clone_genotype,args.output)
create_assignment_file(clone_cells,args.noOfCells,args.op_ass)
