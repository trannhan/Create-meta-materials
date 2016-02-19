#!/bin/bash
#$ -N test6a
#$ -pe mpi-fill 16
#$ -l h_rt=10:00:00
#$ -M nhantran@ksu.edu -m bae
#$ -l mem=500M
# -l ib=TRUE
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
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.000001 -d 0.015 -n 1000000 -p 125 -c 8000 >${JOB_NAME}/output2 2>${JOB_NAME}/error2
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.000001 -d 0.02 -n 1000000 -p 125 -c 8000 >${JOB_NAME}/output3 2>${JOB_NAME}/error3
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.000001 -d 0.025 -n 1000000 -p 125 -c 8000 >${JOB_NAME}/output4 2>${JOB_NAME}/error4
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.000001 -d 0.0095 -n 1000000 -p 125 -c 8000 >${JOB_NAME}/output5 2>${JOB_NAME}/error5

# move the output data to subdirectory
# mv sol.dat ${JOB_NAME}
