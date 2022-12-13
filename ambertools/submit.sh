#!/bin/bash  -l
#SBATCH -J test_set_energies          # Job name
#SBATCH -o o.%j                    # Name of stdout output file
#SBATCH -e e.%j                    # Name of stderr error file
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 08:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=lobo@ucsb.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

cd $SLURM_SUBMIT_DIR
/bin/hostname

python compute_amber_energy.py | tee test_set_output.txt


