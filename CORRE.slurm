#!/bin/bash

#SBATCH --partition 95g
#SBATCH --time 8-00:00:00  ## TOTAL running walltime: min,  min:sec,  h:min:sec, d-h:min:sec
#SBATCH --nodes 1
#SBATCH --ntasks 1
##SBATCH --mem-per-cpu 2000 ## in MB. Alternatively, --mem asks for mem (in MB) *per node*
#SBATCH --mem 80000
#SBATCH --job-name AF_15atom
#SBATCH --mail-type FAIL
#SBATCH --mail-type END
#SBATCH --mail-user crislemon93@gmail.com

export MKL_NUM_THREADS=1

# User variables
SCRATCH="/scratch/${USER}/${SLURM_JOB_ID}"
RESULT_DIR=${SLURM_SUBMIT_DIR}/results-${SLURM_JOB_ID}

# Create result dir
mkdir ${RESULT_DIR}

# Create /scratch directory
mkdir -p ${SCRATCH}

# Copy files to /scratch
cp ${SLURM_SUBMIT_DIR}/{*.py,*.slurm} ${SCRATCH}

# go to scratch
cd ${SCRATCH}

ulimit -s unlimited

########################
# Load environment
########################
module load anaconda/4.3.1-python3.6
python 2D_matrix.py
#python test_AF_FM.py > result-${SLURM_JOB_ID}.out

# cp ${SCRATCH}/* ${SLURM_SUBMIT_DIR}
#cp result-${SLURM_JOB_ID}.out ${SLURM_SUBMIT_DIR}
cp ${SCRATCH}/{*.py,*.slurm,*.txt,*.pdf} ${RESULT_DIR}


