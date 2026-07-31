# Approx. 312 MB for functions etc., 5 MB per iteration; goal: maybe 100-200 iterations, across ~8 cores that's
# 25 iterations/core; ~512 MB should suffice, push to 768 for overhead. Approx. 27 seconds/iteration, for 25 iterations
# 10 minutes required; let's see about expanding the test.

#!/bin/bash
#SBATCH --account=
#SBATCH --time=00:30:00
#SBATCH --mem=512
#SBATCH --job-name=sim_arrays
#SBATCH --output=sim_arrays-%J.out

Rscript sim_arrays_job.R 25 %J