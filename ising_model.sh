#!/bin/bash
#$ -S /bin/bash
#$ -q single.q
#$ -v OMP_NUM_THREADS=1
#$ -cwd
#$ -j y

echo "nodes list"
cat ${PE_HOSTFILE}
pwd
echo "--------------------------------"
echo "jobid:       ${PE_JOBID}"
echo "nodes:       ${NNODES}"
echo "cores:       ${NCORES}"
echo "Nodefile:    ${PE_HOSTFILE}"
echo "# time-start "`date`
echo "--------------------------------"

##source ~/.bashrc

./ising_model
