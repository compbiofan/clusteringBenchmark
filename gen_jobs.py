import argparse 

# Read each line after # in args.input and include them in the slurm file.
# And include the header
def create_job_file(fname,n,mpc,p,cluster_algo):
    with open(fname,'r') as f:
        lines = f.readlines()
    count = 0
        
    op_fname = (fname.split("."))[0]
    for line in lines:
        if "#" in line:
            continue
        line_items = line.split(" ")
        #print(" Line items ",line_items)
        for item in line_items:
            if 'output' in item and cluster_algo == 'scclone':
                opdir_list = (item.strip()).split('/')
                opdir_name = opdir_list[0]+"/"+opdir_list[1]+"/"+opdir_list[2]
            elif 'output' in item:
                #print(" Output ",item)
                opdir_name = item.strip()

        #opdir_name = (line_items[8].strip()).replace(".","")
        #print(opdir_name)
        with open(op_fname+"."+str(count)+".slurm",'w') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --job-name="+cluster_algo+str(count)+"\n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -n "+str(n)+"\n")
            f.write("#SBATCH --mem-per-cpu="+mpc+"\n")
            f.write("#SBATCH -p "+p+"\n")
            f.write("#SBATCH --mail-type=ALL\n")
            f.write("#SBATCH --output="+opdir_name+"/output"+str(count)+".out\n")
            f.write("#SBATCH --error="+opdir_name+"/error"+str(count)+".out\n")
            f.write(line)
        count=count+1

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input .sh file")
parser.add_argument("-n","--n",dest="n", help="No of cpus")
parser.add_argument("-mem_per_cpu", "--mem_per_cpu",dest ="mem_per_cpu", help="Memory per CPU")
parser.add_argument("-p", "--p",dest ="p", help="partition")
parser.add_argument("-algo", "--algo",dest ="algo", help="Which cluster algo you are running")
args = parser.parse_args()

create_job_file(args.input, args.n, args.mem_per_cpu, args.p, args.algo)
