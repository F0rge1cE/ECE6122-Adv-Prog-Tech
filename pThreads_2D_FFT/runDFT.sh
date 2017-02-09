#!/bin/tcsh
#PBS -q class
#PBS -l nodes=1:sixcore
#PBS -l walltime=00:10:00
#PBS -N Xueyang Xu
# The below chnages the working directory to the location of
# your threaded DFT program
cd ThreadsTransform2D
./threadDFT2d



