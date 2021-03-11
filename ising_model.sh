#!/bin/bash
#$ -S /bin/bash
#$ -q single.q
#$ -v OMP_NUM_THREADS=1
#$ -cwd
#$ -j y

pwd
echo "# time-start "`date`
echo "--------------------------------------------"
##source ~/.bashrc
./ising_model
echo "--------------------------------------------"
