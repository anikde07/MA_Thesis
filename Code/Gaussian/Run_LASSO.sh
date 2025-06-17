#!/bin/bash
#SBATCH -J "LASSO"          						# Job name
#SBATCH -o log/Lasso_Polynomial_%A_%a.out       # Output file
#SBATCH -n 1                       				# Number of tasks
#SBATCH -c 4                      				# Number of cores per task
#SBATCH -t 5-21:00:00              				# Time limit
#SBATCH --array=1-9               				# Create jobs in the array

export LC_MONETARY=en_US.utf8 LC_PAPER=en_US.utf8 LC_MEASUREMENT=en_US.utf8
module load linux-rocky8-x86_64/gcc/12.2.0/r/4.4.1-inmpzzy

N_values=(100 300 500 100 300 500 100 300 500)
T_values=(100 100 100 300 300 300 500 500 500)
N=${N_values[$SLURM_ARRAY_TASK_ID-1]}
T=${T_values[$SLURM_ARRAY_TASK_ID-1]}

Rscript LASSO.R $N $T
