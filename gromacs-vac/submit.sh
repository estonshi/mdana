#!/bin/bash

### Job Name
#PBS -N Rec_02
### OutPut Files
#PBS -o Rec_02.stdout
#PBS -e Rec_02.stderr
### Queue Name
#PBS -q low
### Number of nodes
#PBS -l nodes=7:ppn=24

#source /home/ycshi/.bashrc
#source /public/software/profile.d/mpi_openmpi-2.0.0.sh
source /public/software/profile.d/mpi_openmpi-2.0.0-gcc-5.4.0.sh

cd $PBS_O_WORKDIR
MDRUN=/public/home/ycshi/Documents/package/gromacs-4.6.7/bin/mdrun_mpi
mpirun -n 168 $MDRUN -s md.tpr -o traj.trr -c out.gro
