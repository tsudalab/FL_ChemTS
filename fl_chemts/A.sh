#!/bin/sh

#$ -S /bin/bash
#$ -cwd
#$ -pe impi 20
#$ -jc pcc-skl
##$ -jc pcc-large

. /fefs/opt/x86_64/Gaussian/envset.sh
source ~/.bashrc

ulimit -s unlimited

python /home/sumita/GaussianRun_0.8/main.py  CT_sample.sdf  > log
#./test_unix.sh
