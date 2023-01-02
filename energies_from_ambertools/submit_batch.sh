#!/bin/bash  -l
#SBATCH -J training_batch_1-40          # Job name
#SBATCH -o o.%j                    # Name of stdout output file
#SBATCH -e e.%j                    # Name of stderr error file
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 40               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=lobo@ucsb.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

cd $SLURM_SUBMIT_DIR
/bin/hostname

srun --ntasks 1 python compute_amber_energy_batches.py 1  200 | tee  train_batch1_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 2  200 | tee  train_batch2_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 3  200 | tee  train_batch3_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 4  200 | tee  train_batch4_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 5  200 | tee  train_batch5_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 6  200 | tee  train_batch6_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 7  200 | tee  train_batch7_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 8  200 | tee  train_batch8_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 9  200 | tee  train_batch9_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 10 200 | tee train_batch10_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 11 200 | tee train_batch11_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 12 200 | tee train_batch12_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 13 200 | tee train_batch13_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 14 200 | tee train_batch14_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 15 200 | tee train_batch15_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 16 200 | tee train_batch16_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 17 200 | tee train_batch17_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 18 200 | tee train_batch18_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 19 200 | tee train_batch19_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 20 200 | tee train_batch20_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 21 200 | tee train_batch21_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 22 200 | tee train_batch22_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 23 200 | tee train_batch23_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 24 200 | tee train_batch24_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 25 200 | tee train_batch25_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 26 200 | tee train_batch26_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 27 200 | tee train_batch27_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 28 200 | tee train_batch28_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 29 200 | tee train_batch29_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 30 200 | tee train_batch30_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 31 200 | tee train_batch31_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 32 200 | tee train_batch32_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 33 200 | tee train_batch33_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 34 200 | tee train_batch34_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 35 200 | tee train_batch35_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 36 200 | tee train_batch36_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 37 200 | tee train_batch37_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 38 200 | tee train_batch38_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 39 200 | tee train_batch39_output.txt &
srun --ntasks 1 python compute_amber_energy_batches.py 40 200 | tee train_batch40_output.txt &
wait



