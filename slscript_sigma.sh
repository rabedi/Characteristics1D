#!/bin/bash
#This file is a submission script to request the ISAAC resources from Slurm
#SBATCH -J 	m050			       #The name of the job
#SBATCH -A ACF-UTK0011              # The project account to be charged
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=24        # cpus per node #SBATCH --cpus-per-task=1 #SBATCH --mem-per-cpu=17000
#SBATCH --partition=campus-sigma          # If not specified then default is "campus"
#SBATCH --time=1-00:00:00             # Wall time (days-hh:mm:ss)
#SBATCH --error=job.e%J	       # The file where run time errors will be dumped
#SBATCH --output=job.o%J	       # The file where the output of the terminal will be dumped
#SBATCH --qos=campus
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rabedi@utk.edu
####------ ACF mpich ------:

for ((i=0; i<=10; i++))
do
	sh script${i}.sh > output_${i}.txt &
done
wait
############ end of PBSscript ##########
