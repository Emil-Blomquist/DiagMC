#!/bin/bash

#SBATCH -J DiagMC               # A single job name for the array
#SBATCH -t 0-10:00              # Maximum execution time (D-HH:MM)
#SBATCH -n 1                    # Number of cores
#SBATCH -e error_file.e         # Standard error
#SBATCH -o output_file.o        # Standard output

echo ./bin/run $*
./bin/run $*