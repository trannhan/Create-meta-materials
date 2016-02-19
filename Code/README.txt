Compile:
-------------------------------------------------------------------------------------------------------------------------------------
.build
or: 
make


Run:
-------------------------------------------------------------------------------------------------------------------------------------
mpirun -n <number_of_processors> ./<binary_file> -n <size> -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol> -log_summary
or:
mpiexec -n <number_of_processors> ./<binary_file> <options>


Help:
-------------------------------------------------------------------------------------------------------------------------------------
scat -h


Submit job:
-------------------------------------------------------------------------------------------------------------------------------------
qsub run.sh


Requirements: 
-------------------------------------------------------------------------------------------------------------------------------------
MPICH, LAPACK, BLAS, PETSC, GCC


Top iterative methods with respect to least error(from top to bottom):
-------------------------------------------------------------------------------------------------------------------------------------
specest
bcgs
cgs
tfqmr 
cr
gmres, lgmres, fgmres  
cg (runs fastest)
minres
symmlq 


Notes
-------------------------------------------------------------------------------------------------------------------------------------
- use #include 'mpif.h' to work with mpif90.gfortran. Do not use "use mpi"
- Configure PETSC: 
./configure --with-scalar-type=complex --with-precision=double --with-pthreadclasses --with-debugging=no --with-fortran-kernels=1
./configure --with-scalar-type=complex --with-precision=single --with-debugging=no 
(--with-scalar-type=real,complex --with-precision=single,double,longdouble,int,matsingle)
(--with-cc=mpicc or --with-mpi-dir [and not --with-cc=gcc])