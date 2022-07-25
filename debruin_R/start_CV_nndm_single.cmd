#!/bin/bash

# set the number of nodes
#SBATCH --nodes=1

# set the number of CPU cores per node
#SBATCH --ntasks-per-node 1

# How much memory is needed (per node)
#SBATCH --mem=15GB

# set a partition
#SBATCH --partition normal

# set max wallclock time
#SBATCH --time=90:00:00

# set name of job
#SBATCH --job-name=nndm_2ndwave

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# set an output file
#SBATCH --output output_nndm_single_.dat

# send mail to this address
#SBATCH --mail-user=jbahlmann@uni-muenster.de

# run the application
module add palma/2020b
module add foss R GDAL
R CMD BATCH --vanilla CV_nndm_single.R
