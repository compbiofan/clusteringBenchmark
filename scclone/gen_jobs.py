import argparse 

# Read each line after # in args.input and include them in the slurm file.
# And include the header
def create_job_file(fname,n,mpc,p):
    with open(fname,'r') as f:
        lines = f.readlines()
   
    count = 0
    op_fname = (fname.split("."))[0]
    for line in lines:
        if "#" in line:
            continue
        line_items = line.split(" ")
        opdir_name = (((line_items[9].strip()).replace(".","")).split("/"))[2]
        #print(opdir_name)
        with open(op_fname+"."+str(count)+".slurm",'w') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --job-name=scclone_job"+str(count)+"\n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -n "+str(n)+"\n")
            f.write("#SBATCH --mem-per-cpu="+mpc+"\n")
            f.write("#SBATCH -p "+p+"\n")
            f.write("#SBATCH --mail-type=ALL\n")
            f.write("#SBATCH --output=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/"+opdir_name+"/scclone_output"+str(count)+".out\n")
            f.write("#SBATCH --error=/gpfs/research/fangroup/rk18g/longitudinal/clustering/scclone/output/"+opdir_name+"/scclone_error"+str(count)+".out\n")
            f.write(line)
        count=count+1

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input .sh file")
parser.add_argument("-n","--n",dest="n", help="No of cpus")
parser.add_argument("-mem_per_cpu", "--mem_per_cpu",dest ="mem_per_cpu", help="Memory per CPU")
parser.add_argument("-p", "--p",dest ="p", help="partition")
args = parser.parse_args()

create_job_file(args.input, args.n, args.mem_per_cpu, args.p)
