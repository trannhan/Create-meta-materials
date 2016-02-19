
/* Program usage:  mpiexec -n <procs> scat [-help] [all PETSc options] */

static char help[] = "Solves particle-scattering problem in parallel.\n\
References:\n\
[595]  A.G. Ramm,  Wave scattering by many small bodies and creating materials with a desired refraction coefficient, Afrika Matematika, 22, N1, (2011), 33-55.\n\
Input parameters include:\n\
  -a <particle_radius>         : particle radius\n\
  -d <particles_distance>      : distance between neighboring particles, 1 >> d >> a, default value cuberoot[a^(2-Kappa)]\n\
  -n <total_particles>         : total number of particles, n = O(1/(a^(2-Kappa)))\n\
  -p <total_cubes>             : total number of embeded small cubes in the domain containing all particles (for solving the reduced system)\n\
  -c <total_collocation_points>: total number of collocation points in the domain containing all particles (for solving the integral equation)\n\
  -k <kappa>                   : power constant with respect to the radius of particles, kappa is in (0,1), default value 0.99\n\
  -vol <volume>                : volume of the domain that contains all particles, default value 1\n\
  -ori <original_refraction>   : original refraction coefficient, default value 1\n\
  -des <desired_refraction>    : desired refraction coefficient, default value sqrt(0.2)\n\
  -dis <distribution>          : distribution of particles, default value 1 for uniform distribution\n\
  -view_RHS                    : write RHS vector to stdout\n\
  -view_solution               : write solution vector to stdout\n\n";

/*
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/

//#include <petscksp.h>
#include "scattering3DS.h"
#include "scattering3DP.h"
#include "scatteringIE_Riemann.h"
#include "scatteringIE_Collocation.h"


using namespace std;

extern inline const PetscErrorCode ScatSMatMult(Mat,Vec,Vec);
extern inline const PetscErrorCode ScatPMatMult(Mat,Vec,Vec);
extern inline const PetscErrorCode ScatIEMatMult(Mat,Vec,Vec);
extern PetscErrorCode DiffP(Vec,Vec,Vec,Vec);
extern PetscErrorCode DiffPCollocation(Vec,Vec,Vec,Vec);
extern PetscErrorCode DiffIE(Vec,Vec,Vec,Vec);


Scattering3DS<std::vector<PetscReal>,std::vector<PetscScalar>,PetscScalar,PetscReal,PetscInt,PetscInt> Scat3DS;
Scattering3DP<std::vector<PetscReal>,std::vector<PetscScalar>,PetscScalar,PetscReal,PetscInt,PetscInt> Scat3DP;
ScatteringIE_Riemann<std::vector<PetscReal>,std::vector<PetscScalar>,PetscScalar,PetscReal,PetscInt,PetscInt> ScatIE;


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x,y,z,r,rP,rIE,rC,u0,u0P,u0IE,Distance,nParticlePerCube;  /* approx solution, residual, RHS, solution-difference */
  Mat            A,B,IE;           		                          /* linear system matrix */
  KSP            ksp,kspP,kspIE;      	                                 /* linear solver context */
  //PC             pc;          		                           	     /* Preconditioner context */
  PetscReal      norm;          		                                 /* norm of solution error */
  PetscInt       its;           		                                 /* number of iterations reached */
  //PetscInt       Istart,Iend;
  PetscErrorCode ierr;
  KSPConvergedReason reason;
  
  //Scattering3D:
  PetscInt       n = 0;                 				   //Total particles
  PetscReal      a = 0;                 				   //Particle radius
  PetscReal      ParticleDistance = 0;
  PetscReal      VolQ = 0;              				   //Volume of the domain Q that contains all particles
  PetscReal      Kappa = 0;             				   //Power const with respect to the radius of particles: Kappa in [0,1]
  std::vector<PetscReal> WaveDirection(3,0);WaveDirection[0] = 1; // WaveDirection is a unit vector that indicates the direction of plane wave
  PetscScalar    OriginalRefractionCoef = 0;
  PetscScalar    DesiredRefractionCoef = 0;
  PetscReal      Distribution = 0;
  PetscInt       TotalCubes = 0;
  PetscInt       N = 0;                 				   //Number of sub domains used for solving IE
  PetscScalar    tmp;

  int            rank;
  int            size;
  time_t         time1, time2; //Time usage 
  PetscBool      flg = PETSC_FALSE;

#if defined(PETSC_USE_LOG)
  PetscLogStage  stage;
#endif    

  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  //printf ("Hello from task %d!\n", rank);

  time(&time1);

  ierr = PetscOptionsGetReal(PETSC_NULL,"-a",&a,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-d",&ParticleDistance,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-p",&TotalCubes,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-c",&N,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-k",&Kappa,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-vol",&VolQ,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-ori",&OriginalRefractionCoef,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-des",&DesiredRefractionCoef,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-dis",&Distribution,PETSC_NULL);CHKERRQ(ierr);

  //Scatering init 
  ValidateInput(Kappa,VolQ,a,ParticleDistance,n,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,TotalCubes,N);
  Scat3DS.Input(a,Kappa,WaveDirection,ParticleDistance,n,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,VolQ);
  Scat3DP.Input(a,Kappa,WaveDirection,ParticleDistance,n,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,VolQ,TotalCubes,N); 
  ScatIE.Input(a,Kappa,WaveDirection,ParticleDistance,n,OriginalRefractionCoef,DesiredRefractionCoef,Distribution,VolQ,N); 
  Scat3DS.Init();
  Scat3DP.Init();
  ScatIE.Init();
  if(rank==0)
  {
      Output<std::vector<PetscReal>,PetscScalar,PetscReal,PetscInt,PetscInt>(Kappa,VolQ,a,ParticleDistance,n,WaveDirection,\
                                                                             OriginalRefractionCoef,DesiredRefractionCoef,Distribution,TotalCubes,Scat3DS.BoundaryImpedance,N);
  }
  
  //----------------------------------------SOLVING THE ORIGINAL SYSTEM S-------------------------------------------------------
  ierr = PetscLogStageRegister("Original system", &stage);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stage);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = u0.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */
  ierr = MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,n,n,(void*)0,&A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))ScatSMatMult);CHKERRQ(ierr);  

  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
  //ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
  //ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);

  /*
     Set up RHS and solution vectors
  */
  //ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,&u0);CHKERRQ(ierr);  //Better create vector thread
  ierr = VecCreate(PETSC_COMM_WORLD,&u0);CHKERRQ(ierr);
  ierr = VecSetSizes(u0,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u0);CHKERRQ(ierr);
  ierr = VecDuplicate(u0,&x);CHKERRQ(ierr);

  //Testing:
  //ierr = VecSet(u0,1.0);CHKERRQ(ierr);
  //ierr = VecSet(x,1.0);CHKERRQ(ierr);
  //MatMult(A,x,u0);

  for(int s = 0; s <n ; s++)
  {
      tmp = Scat3DS.InitField(s);
      VecSetValues(u0,1,&s,&tmp,INSERT_VALUES);
  }
  ierr = VecAssemblyBegin(u0);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u0);CHKERRQ(ierr);

  /*
     View the RHS vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_RHS",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVector Right-hand-side of the original linear system:\n");CHKERRQ(ierr);
    ierr = VecView(u0,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */
  ierr = KSPSetType(ksp,KSPSPECEST);CHKERRQ(ierr);
  //ierr = KSPGMRESSetRestart(ksp,10);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1e-3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  //ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  //ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  //ierr = KSPSetUp(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = KSPSolve(ksp,u0,x);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Check the error
  */
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
  ierr = MatMult(A,x,r); CHKERRQ(ierr);
  ierr = VecAXPY(r,-1.0,u0);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  /* Scale the norm */
  /*  norm *= sqrt(1.0/n); */

  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  ierr = KSPGetConvergedReason(ksp,&reason);CHKERRQ(ierr);
  if (reason<0) {
     ierr = PetscPrintf(PETSC_COMM_WORLD,"The original linear system is divergent!\n");CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSolving the original linear system:\nNorm of error:\t%G\nIterations:\t%D\n",norm,its);CHKERRQ(ierr);

  /*
     View the Solution vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_solution",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVector Solution of the original linear system:\n");CHKERRQ(ierr);
    ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  ierr = PetscLogStagePop();CHKERRQ(ierr);

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */  
  ierr = VecDestroy(&u0);CHKERRQ(ierr);
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  //----------------------------------------SOLVING THE REDUCED SYSTEM P-------------------------------------------------------
  ierr = PetscLogStageRegister("Reduced system", &stage);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stage);CHKERRQ(ierr);

  ierr = MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,TotalCubes,TotalCubes,(void*)0,&B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);
  ierr = MatShellSetOperation(B,MATOP_MULT,(void(*)(void))ScatPMatMult);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&u0P);CHKERRQ(ierr);
  ierr = VecSetSizes(u0P,PETSC_DECIDE,TotalCubes);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u0P);CHKERRQ(ierr);
  ierr = VecDuplicate(u0P,&y);CHKERRQ(ierr);

  for(int s = 0; s <TotalCubes ; s++)
  {
      tmp = Scat3DP.InitField(s);
      VecSetValues(u0P,1,&s,&tmp,INSERT_VALUES);
  }
  ierr = VecAssemblyBegin(u0P);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u0P);CHKERRQ(ierr);

  /*
     View the RHS vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_RHS",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVector Right-hand-side of the reduced linear system:\n");CHKERRQ(ierr);
    ierr = VecView(u0P,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
  }

  ierr = KSPCreate(PETSC_COMM_WORLD,&kspP);CHKERRQ(ierr);
  ierr = KSPSetOperators(kspP,B,B,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetType(kspP,KSPSPECEST);CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspP,1e-3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspP);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPSolve(kspP,u0P,y);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Check the error
  */
  ierr = VecDuplicate(y,&rP);CHKERRQ(ierr);
  ierr = MatMult(B,y,rP); CHKERRQ(ierr);
  ierr = VecAXPY(rP,-1.0,u0P);CHKERRQ(ierr);
  ierr = VecNorm(rP,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(kspP,&its);CHKERRQ(ierr);
  /* Scale the norm */
  /*  norm *= sqrt(1.0/TotalCubes); */

  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  ierr = KSPGetConvergedReason(kspP,&reason);CHKERRQ(ierr);
  if (reason<0) {
     ierr = PetscPrintf(PETSC_COMM_WORLD,"The reduced linear system is divergent!\n");CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSolving the reduced linear system:\nNorm of error:\t%G\nIterations:\t%D\n",norm,its);CHKERRQ(ierr);

  /*
     View the Solution vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_solution",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVector Solution of the reduced linear system:\n");CHKERRQ(ierr);
    ierr = VecView(y,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  ierr = PetscLogStagePop();CHKERRQ(ierr); 

  ierr = VecDestroy(&rP);CHKERRQ(ierr);
  ierr = VecDestroy(&u0P);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = KSPDestroy(&kspP);CHKERRQ(ierr);

//----------------------------------------SOLVING THE INTEGRAL EQUATION-------------------------------------------------------

  ierr = PetscLogStageRegister("Integral Equation", &stage);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stage);CHKERRQ(ierr);

  ierr = MatCreateShell(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,(void*)0,&IE);CHKERRQ(ierr);
  ierr = MatSetUp(IE);CHKERRQ(ierr);
  ierr = MatShellSetOperation(IE,MATOP_MULT,(void(*)(void))ScatIEMatMult);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&u0IE);CHKERRQ(ierr);
  ierr = VecSetSizes(u0IE,PETSC_DECIDE,N);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u0IE);CHKERRQ(ierr);
  ierr = VecDuplicate(u0IE,&z);CHKERRQ(ierr);

  for(int s = 0; s <N ; s++)
  {
      tmp = ScatIE.InitField(s);
      VecSetValues(u0IE,1,&s,&tmp,INSERT_VALUES);
  }
  ierr = VecAssemblyBegin(u0IE);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(u0IE);CHKERRQ(ierr);

  /*
     View the RHS vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_RHS",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVector Right-hand-side of the integral equation:\n");CHKERRQ(ierr);
    ierr = VecView(u0IE,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    
  }

  ierr = KSPCreate(PETSC_COMM_WORLD,&kspIE);CHKERRQ(ierr);
  ierr = KSPSetOperators(kspIE,IE,IE,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetType(kspIE,KSPSPECEST);CHKERRQ(ierr);
  ierr = KSPSetTolerances(kspIE,1e-3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(kspIE);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = KSPSolve(kspIE,u0IE,z);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Check the error
  */
  ierr = VecDuplicate(z,&rIE);CHKERRQ(ierr);
  ierr = MatMult(IE,z,rIE); CHKERRQ(ierr);
  ierr = VecAXPY(rIE,-1.0,u0IE);CHKERRQ(ierr);
  ierr = VecNorm(rIE,NORM_2,&norm);CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(kspIE,&its);CHKERRQ(ierr);
  /* Scale the norm */
  /*  norm *= sqrt(1.0/N); */

  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  ierr = KSPGetConvergedReason(kspIE,&reason);CHKERRQ(ierr);
  if (reason<0) {
     ierr = PetscPrintf(PETSC_COMM_WORLD,"The integral equation has no solution!\n");CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nSolving the integral equation:\nNorm of error:\t%G\nIterations:\t%D\n",norm,its);CHKERRQ(ierr);

  /*
     View the Solution vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_solution",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nVector Solution of the integral equation:\n");CHKERRQ(ierr);
    ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  ierr = PetscLogStagePop();CHKERRQ(ierr); 

  ierr = VecDestroy(&rIE);CHKERRQ(ierr);
  ierr = VecDestroy(&u0IE);CHKERRQ(ierr);
  ierr = MatDestroy(&IE);CHKERRQ(ierr);
  ierr = KSPDestroy(&kspIE);CHKERRQ(ierr);


  //----------------------------------------COMPARE SOLUTIONS BETWEEN S, P SYSTEMS WITH IE-------------------------------------------------------

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nCompare solutions:");CHKERRQ(ierr);
  //x: solution to S
  //y: solution to P, size(y) << size(x)
  //z: solution to IE, size(y) <= size(z) <= size(x)

  //Compare S & P:
  ierr = VecDuplicate(y,&Distance);CHKERRQ(ierr);
  ierr = VecDuplicate(y,&nParticlePerCube);CHKERRQ(ierr);
  ierr = VecDuplicate(y,&rC);CHKERRQ(ierr);

  ierr = DiffP(x,y,Distance,nParticlePerCube);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(rC,Distance,nParticlePerCube);CHKERRQ(ierr);
  ierr = VecNorm(rC,NORM_INFINITY,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Original system)   vs (Reduced system):\t%G",norm);CHKERRQ(ierr);

  ierr = VecDestroy(&Distance);CHKERRQ(ierr);
  ierr = VecDestroy(&nParticlePerCube);CHKERRQ(ierr);
  ierr = VecDestroy(&rC);CHKERRQ(ierr);

  //Compare S & IE:
  ierr = VecDuplicate(z,&Distance);CHKERRQ(ierr);
  ierr = VecDuplicate(z,&nParticlePerCube);CHKERRQ(ierr);
  ierr = VecDuplicate(z,&rC);CHKERRQ(ierr);

  ierr = DiffIE(x,z,Distance,nParticlePerCube);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(rC,Distance,nParticlePerCube);CHKERRQ(ierr);
  ierr = VecNorm(rC,NORM_INFINITY,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Integral equation) vs (Original system):\t%G",norm);CHKERRQ(ierr);

  ierr = VecDestroy(&Distance);CHKERRQ(ierr);
  ierr = VecDestroy(&nParticlePerCube);CHKERRQ(ierr);
  ierr = VecDestroy(&rC);CHKERRQ(ierr);

  //Compare P & IE:
  ierr = VecDuplicate(y,&Distance);CHKERRQ(ierr);
  ierr = VecDuplicate(y,&nParticlePerCube);CHKERRQ(ierr);
  ierr = VecDuplicate(y,&rC);CHKERRQ(ierr);

  ierr = DiffPCollocation(z,y,Distance,nParticlePerCube);CHKERRQ(ierr);
  ierr = VecPointwiseDivide(rC,Distance,nParticlePerCube);CHKERRQ(ierr);
  ierr = VecNorm(rC,NORM_INFINITY,&norm);CHKERRQ(ierr);
  //VecAXPY(z,-1,y);
  //ierr = VecNorm(z,NORM_INFINITY,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n(Integral equation) vs (Reduced system):\t%G\n",norm);CHKERRQ(ierr);

  ierr = VecDestroy(&Distance);CHKERRQ(ierr);
  ierr = VecDestroy(&nParticlePerCube);CHKERRQ(ierr);
  ierr = VecDestroy(&rC);CHKERRQ(ierr);

  //Check total time used:
  time(&time2); 
  if(rank==0)
  {
      checkTime(time1,time2);
  }

  ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&y);CHKERRQ(ierr);
  ierr = VecDestroy(&z);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  PetscFinalize();
  
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
/*
#undef __FUNCT__
#define __FUNCT__ "MatMult"
//scatter the whole vector xx
PetscErrorCode MatMult(Mat mat,Vec xx,Vec yy)
{
  PetscInt       n,s,t,ystart,yend;
  PetscScalar    tmp,v;
  PetscErrorCode ierr;
  PetscScalar    Green=2,h=1,a2k=1;
  Vec            vec;
  VecScatter     scat;
  IS             is;

  PetscFunctionBegin;
  ierr  =  VecGetOwnershipRange(yy,&ystart,&yend);CHKERRQ(ierr);

  ierr = VecGetSize(xx,&n);CHKERRQ(ierr);

  VecCreate(PETSC_COMM_SELF,&vec);
  VecSetSizes(vec,PETSC_DECIDE,n);
  VecSetFromOptions(vec);

  ISCreateStride(PETSC_COMM_WORLD,n,0,1,&is);
  VecScatterCreate(xx,PETSC_NULL,vec,is,&scat);
  //VecScatterCreateToAll(xx,&scat,&vec);
  VecScatterBegin(scat,xx,vec,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(scat,xx,vec,INSERT_VALUES,SCATTER_FORWARD);

  for(s = ystart; s <yend ; s++)
  {
      tmp = 0;
      for(t = 0; t < n; t++)
      {
          VecGetValues(vec,1,&t,&v);
          tmp += Scat3DS.CoefMat(s,t)*v;
      }
      VecSetValues(yy,1,&s,&tmp,INSERT_VALUES);
  }

  VecDestroy(&vec);
  VecScatterDestroy(&scat);
  ISDestroy(&is);

  PetscFunctionReturn(0);
}
*/

////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "ScatSMatMult"
inline const PetscErrorCode ScatSMatMult(Mat mat,Vec xx,Vec yy)
{
  PetscInt       n,s,t;
  PetscScalar    tmp,v;
  //PetscErrorCode ierr;
  PetscInt       xstart,xend;
  const char     *prefix;

  PetscFunctionBegin;
  VecGetOwnershipRange(xx,&xstart,&xend);
  VecGetSize(yy,&n);

  VecGetOptionsPrefix(yy,&prefix);
  if(prefix!="ZeroSetDone")
  {
      tmp = 0;
      for(s = 0; s <n ; s++)
      {
         VecSetValues(yy,1,&s,&tmp,INSERT_VALUES);
      }
      VecAssemblyBegin(yy);
      VecAssemblyEnd(yy);
      VecSetOptionsPrefix(yy,"ZeroSetDone");
  }

  for(s = 0; s <n ; s++)
  {
      tmp = 0;
      for(t = xstart; t < xend; t++)
      {
          VecGetValues(xx,1,&t,&v);
          tmp += Scat3DS.CoefMatFast(s,t)*v;
      }
      VecSetValues(yy,1,&s,&tmp,ADD_VALUES);
  }
  VecAssemblyBegin(yy);
  VecAssemblyEnd(yy);

  PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "ScatPMatMult"
inline const PetscErrorCode ScatPMatMult(Mat mat,Vec xx,Vec yy)
{
  PetscInt       n,s,t;
  PetscScalar    tmp,v;
  //PetscErrorCode ierr;
  PetscInt       xstart,xend;
  const char     *prefix;

  PetscFunctionBegin;
  VecGetOwnershipRange(xx,&xstart,&xend);
  VecGetSize(yy,&n);

  VecGetOptionsPrefix(yy,&prefix);
  if(prefix!="ZeroSetDone")
  {
      tmp = 0;
      for(s = 0; s <n ; s++)
      {
         VecSetValues(yy,1,&s,&tmp,INSERT_VALUES);
      }
      VecAssemblyBegin(yy);
      VecAssemblyEnd(yy);
      VecSetOptionsPrefix(yy,"ZeroSetDone");
  }

  for(s = 0; s <n ; s++)
  {
      tmp = 0;
      for(t = xstart; t < xend; t++)
      {
          VecGetValues(xx,1,&t,&v);
          tmp += Scat3DP.CoefMatFast(s,t)*v;
      }
      VecSetValues(yy,1,&s,&tmp,ADD_VALUES);
  }
  VecAssemblyBegin(yy);
  VecAssemblyEnd(yy);

  PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "ScatIEMatMult"
inline const PetscErrorCode ScatIEMatMult(Mat mat,Vec xx,Vec yy)
{
  PetscInt       n,s,t;
  PetscScalar    tmp,v;
  //PetscErrorCode ierr;
  PetscInt       xstart,xend;
  const char     *prefix;

  PetscFunctionBegin;
  VecGetOwnershipRange(xx,&xstart,&xend);
  VecGetSize(yy,&n);

  VecGetOptionsPrefix(yy,&prefix);
  if(prefix!="ZeroSetDone")
  {
      tmp = 0;
      for(s = 0; s <n ; s++)
      {
         VecSetValues(yy,1,&s,&tmp,INSERT_VALUES);
      }
      VecAssemblyBegin(yy);
      VecAssemblyEnd(yy);
      VecSetOptionsPrefix(yy,"ZeroSetDone");
  }

  for(s = 0; s <n ; s++)
  {
      tmp = 0;
      for(t = xstart; t < xend; t++)
      {
          VecGetValues(xx,1,&t,&v);
          tmp += ScatIE.CoefMat(s,t)*v;
      }
      VecSetValues(yy,1,&s,&tmp,ADD_VALUES);
  }
  VecAssemblyBegin(yy);
  VecAssemblyEnd(yy);

  PetscFunctionReturn(0);
}


////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "DiffP"
PetscErrorCode DiffP(Vec originalVec, Vec reducedVec, Vec Distance, Vec nParticlePerCube)
{
    // Find the distance between 2 solutions originalVec and reducedVec
    PetscInt    s,t;
    PetscInt    n,xstart,xend;
    PetscScalar tmp,v,diff;
    Vec         vec;
    VecScatter  scat;
    IS          is;

    PetscFunctionBegin;
    VecGetOwnershipRange(originalVec,&xstart,&xend);
    VecGetSize(reducedVec,&n);

    VecCreate(PETSC_COMM_SELF,&vec);
    VecSetSizes(vec,PETSC_DECIDE,n);
    VecSetFromOptions(vec);

    ISCreateStride(PETSC_COMM_WORLD,n,0,1,&is);
    VecScatterCreate(reducedVec,PETSC_NULL,vec,is,&scat);
    //VecScatterCreateToAll(reducedVec,&scat,&vec);
    VecScatterBegin(scat,reducedVec,vec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scat,reducedVec,vec,INSERT_VALUES,SCATTER_FORWARD);

    for(s=xstart;s<xend;s++)
    {
        t = Scat3DP.FindCube(s);  
        VecGetValues(vec,1,&t,&tmp);
        VecGetValues(originalVec,1,&s,&v);

        diff = v - tmp;
        VecSetValues(Distance,1,&t,&diff,ADD_VALUES);
        tmp = 1;
        VecSetValues(nParticlePerCube,1,&t,&tmp,ADD_VALUES);        
    }   

    VecDestroy(&vec);
    VecScatterDestroy(&scat);
    ISDestroy(&is);

    PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "DiffPCollocation"
PetscErrorCode DiffPCollocation(Vec originalVec, Vec reducedVec, Vec Distance, Vec nCollocationPointsPerCube)
{
    // Find the distance between 2 solutions originalVec and reducedVec
    PetscInt    s,t;
    PetscInt    n,xstart,xend;
    PetscScalar tmp,v,diff;
    Vec         vec;
    VecScatter  scat;
    IS          is;

    PetscFunctionBegin;
    VecGetOwnershipRange(originalVec,&xstart,&xend);
    VecGetSize(reducedVec,&n);

    VecCreate(PETSC_COMM_SELF,&vec);
    VecSetSizes(vec,PETSC_DECIDE,n);
    VecSetFromOptions(vec);

    ISCreateStride(PETSC_COMM_WORLD,n,0,1,&is);
    VecScatterCreate(reducedVec,PETSC_NULL,vec,is,&scat);
    //VecScatterCreateToAll(reducedVec,&scat,&vec);
    VecScatterBegin(scat,reducedVec,vec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scat,reducedVec,vec,INSERT_VALUES,SCATTER_FORWARD);

    for(s=xstart;s<xend;s++)
    {
        t = Scat3DP.FindCubeOfCollocation(s);  
        VecGetValues(vec,1,&t,&tmp);
        VecGetValues(originalVec,1,&s,&v);

        diff = v - tmp;
        VecSetValues(Distance,1,&t,&diff,ADD_VALUES);
        tmp = 1;
        VecSetValues(nCollocationPointsPerCube,1,&t,&tmp,ADD_VALUES);        
    }   

    VecDestroy(&vec);
    VecScatterDestroy(&scat);
    ISDestroy(&is);

    PetscFunctionReturn(0);
}

////////////////////////////////////////////////////////////////////////////////////

#undef __FUNCT__
#define __FUNCT__ "DiffIE"
PetscErrorCode DiffIE(Vec originalVec, Vec reducedVec, Vec Distance, Vec nParticlePerCube)
{
    // Find the distance between 2 solutions originalVec and reducedVec
    PetscInt    s,t;
    PetscInt    n,xstart,xend;
    PetscScalar tmp,v,diff;
    Vec         vec;
    VecScatter  scat;
    IS          is;

    PetscFunctionBegin;
    VecGetOwnershipRange(originalVec,&xstart,&xend);
    VecGetSize(reducedVec,&n);

    VecCreate(PETSC_COMM_SELF,&vec);
    VecSetSizes(vec,PETSC_DECIDE,n);
    VecSetFromOptions(vec);

    ISCreateStride(PETSC_COMM_WORLD,n,0,1,&is);
    VecScatterCreate(reducedVec,PETSC_NULL,vec,is,&scat);
    //VecScatterCreateToAll(reducedVec,&scat,&vec);
    VecScatterBegin(scat,reducedVec,vec,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scat,reducedVec,vec,INSERT_VALUES,SCATTER_FORWARD);

    for(s=xstart;s<xend;s++)
    {
        t = ScatIE.FindCube(s);  
        VecGetValues(vec,1,&t,&tmp);
        VecGetValues(originalVec,1,&s,&v);

        diff = v - tmp;
        VecSetValues(Distance,1,&t,&diff,ADD_VALUES);
        tmp = 1;
        VecSetValues(nParticlePerCube,1,&t,&tmp,ADD_VALUES);        
    }   
  
    VecDestroy(&vec);
    VecScatterDestroy(&scat);
    ISDestroy(&is);

    PetscFunctionReturn(0);
}

