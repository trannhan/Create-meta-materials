#!/bin/bash
#$ -N norm
#$ -pe mpi-fill 16
#$ -l h_rt=4:00:00
#$ -M nhantran@ksu.edu -m bae
#$ -l mem=1G
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
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.0001 -d 0.05 -n 10000 -p 125 -c 8000 -ksp_compute_eigenvalues -pc_type none -ksp_monitor_singular_value -ksp_gmres_restart 1000 >${JOB_NAME}/output0 2>${JOB_NAME}/error0
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.00001 -d 0.023 -n 100000 -p 125 -c 8000 -ksp_compute_eigenvalues -pc_type none -ksp_monitor_singular_value -ksp_gmres_restart 1000 >${JOB_NAME}/output1 2>${JOB_NAME}/error1
mpirun -np $NSLOTS ./scat -ksp_type gmres  -a 0.000001 -d 0.01 -n 1000000 -p 125 -c 8000 -ksp_compute_eigenvalues -pc_type none -ksp_monitor_singular_value -ksp_gmres_restart 1000 >${JOB_NAME}/output2 2>${JOB_NAME}/error2

# move the output data to subdirectory
# mv sol.dat ${JOB_NAME}
