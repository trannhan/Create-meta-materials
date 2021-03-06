#!/bin/bash
#$ -N scat3
#$ -pe mpi-fill 8
#$ -l h_rt=1:00:00
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
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.1 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output1 2>${JOB_NAME}/error1
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.15 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output2 2>${JOB_NAME}/error2
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.2 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output3 2>${JOB_NAME}/error3
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.25 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output4 2>${JOB_NAME}/error4
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.3 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output5 2>${JOB_NAME}/error5
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.35 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output6 2>${JOB_NAME}/error6
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.4 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output7 2>${JOB_NAME}/error7
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.45 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output8 2>${JOB_NAME}/error8
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.08 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output9 2>${JOB_NAME}/error9
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.06 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output10 2>${JOB_NAME}/error10
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.04 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output11 2>${JOB_NAME}/error11
#mpirun -np $NSLOTS ./scat -a 0.001 -d 0.02 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output12 2>${JOB_NAME}/error12
mpirun -np $NSLOTS ./scat -a 0.001 -d 0.005 -n 1000 -p 125 -c 1000 >${JOB_NAME}/output13 2>${JOB_NAME}/error13

# move the output data to subdirectory
# mv sol.dat ${JOB_NAME}
