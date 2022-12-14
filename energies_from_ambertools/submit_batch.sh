#!/bin/bash  -l
#SBATCH -J training_batch1_energies          # Job name
#SBATCH -o o.%j                    # Name of stdout output file
#SBATCH -e e.%j                    # Name of stderr error file
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 18:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=lobo@ucsb.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

cd $SLURM_SUBMIT_DIR
/bin/hostname

python compute_amber_energy_batches.py 1 | tee train_batch1_output.txt


