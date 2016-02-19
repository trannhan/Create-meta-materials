#!/bin/bash
#$ -N test5a
#$ -pe mpi-fill 16
#$ -l h_rt=01:00:00
#$ -M nhantran@ksu.edu -m bae
#$ -l mem=200M
#$ -cwd
# -o output


###############################################################################
# Notes for Sun Grid Engine Lines
###############################################################################
# 
# 1) Change the line '-N ...' to a reasonable job name.
# 2) Change the line '-pe ...' to use a good number of CPUs.
# 3) Change the line '-l h_rt ..." to a good overestimate of the runtime.
# 4) Change the line '-M ...' to use the correct email address.
#
###############################################################################


# change directory to working directory
cd ${SGE_O_WORKDIR}

# make a subdirectory to place data and output in
mkdir -p ${JOB_NAME}

# copy the job script to the subdirectory
cp ${JOB_NAME}.sh ${JOB_NAME}

# run the simulation (NOTE: change parameters here)
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.009 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output7 2>${JOB_NAME}/error7
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.011 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output8 2>${JOB_NAME}/error8
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.012 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output9 2>${JOB_NAME}/error9
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.013 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output10 2>${JOB_NAME}/error10
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.014 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output11 2>${JOB_NAME}/error11
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.016 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output12 2>${JOB_NAME}/error12
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.02 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output13 2>${JOB_NAME}/error13
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.025 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output14 2>${JOB_NAME}/error14
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.03 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output15 2>${JOB_NAME}/error15
#mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.035 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output16 2>${JOB_NAME}/error16
mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.023 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output17 2>${JOB_NAME}/error17
mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.024 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output18 2>${JOB_NAME}/error18
mpirun -np $NSLOTS ./scat -ksp_type gmres -a 0.00001 -d 0.026 -n 100000 -p 125 -c 8000 >${JOB_NAME}/output19 2>${JOB_NAME}/error19

# move the output data to subdirectory
# mv sol.dat ${JOB_NAME}
