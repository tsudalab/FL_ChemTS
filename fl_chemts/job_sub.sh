#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -pe impi 30 #minimum 2
##$ -par 20 -pe impi 1200
##$ -mods l_hard h_vmem 9G
##$ -mods l_hard mem_req 9G
#$ -jc pcc-skl.168h #pcc-normal.72h # pcc-skl.168h
## Load the Intel MPI environment variables
. /fefs/opt/x86_64/Gaussian/envset.sh

ulimit -s unlimited

#g16 H2O
#. /fefs/opt/x86_64/intel/parallel_studio_xe_2017/impi/2017.2.174/bin64/mpivars.sh

source /home/terayama/.bashrc

#export PATH="/home/yang/anaconda2/bin:$PATH"
#export PATH="/home/terayama/anaconda3/bin:$PATH"

source activate py2
export KERAS_BACKEND=tensorflow
## -pe impi argument 2161 is automatically transferred to the number of executable processes. Note, however, that -bootstrap sge is required.
mpiexec -bootstrap sge -n 30 python mpi_thread_chemts_tree_vl_0218.py
source deactivate
