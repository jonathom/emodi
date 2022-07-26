#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node 10

# How much memory is needed (per node)
#SBATCH --mem=100GB

# set a partition
#SBATCH --partition normal

# set max wallclock time
#SBATCH --time=20:00:00

# set name of job
#SBATCH --job-name=random_emodi

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# set an output file
#SBATCH --output output_random.dat

# send mail to this address
#SBATCH --mail-user=jbahlmann@uni-muenster.de

# run the application
module add palma/2020b
module add foss R GDAL
R CMD BATCH --vanilla CV_random_flds.R
