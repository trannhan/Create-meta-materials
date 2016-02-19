#ifndef COMMON_H
#define COMMON_H

#include <petscksp.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <complex.h>
#include <stdio.h>
//#include <string>
//#include <sstream>


#define _USE_MATH_DEFINES
#define ROOT 0
#define ZERO 1.0e-6
#define PI   (3.14159265359)
#define PI4  (M_PI*4)


//In Optics
//#define C   (3*1.0e+10)   // Speed of light cm/s
//#define F   (1.0e+14)     // Frequency Hz
//#define K   (2*PI*F/C)    // Wave number K = 2PI/lambda

//Acoustic waves
#define C   (34400)       // Speed of light cm/s
#define F   (1000)        // Frequency Hz
#define K   (2*PI*F/C)    // Wave number K = 2PI/lambda


using namespace std;


///////////////////////////////////////////////////////////////////

template<class Vector,class Real>
Real Norm(Vector v, int norm)
{
    Real sum = 0;
    int n = v.size();

    for(int i=0;i<n;i++)
    {
        sum += pow(v[i],norm);
    }
    return pow(sum,1.0/norm);
}

template<class Real, class Vector>
Real Dot(Vector v, Vector u)
{
    Real sum = 0;
    int n = v.size();

    for(int i=0;i<n;i++)
    {
        sum += v[i]*u[i];
    }
    return sum;
}

template<class Real, class Complex>
Real cnorm(Complex c)
{
    return sqrt(creal(c)*creal(c)+cimag(c)*cimag(c));
}

bool is_nth_power(int a, int n) 
{
  if(n <= 0)
    return false;
  if((a < 0) && (n % 2 == 0))
    return false;
  a = abs(a);

  int b = pow(a, 1. / n);
  return (pow((double) b, n) == a || pow((double) (b+1), n) == a);
}

template <class Complex, class Real, class Long, class Integer> 
void ValidateInput(Real& Kappa, Real& VolQ, Real& ParRadius, Real& ParticleDistance, Long& TotalParticles,\
                  Complex& OriginalRefractionCoef, Complex& DesiredRefractionCoef, Real& Distribution, Integer& TotalCubes, Integer& TotalCollocationPoints)
{
  if( (Kappa<=0) || (Kappa>=1) )
  {
     Kappa = 0.99;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Kappa should be in (0,1). Set the default value!");
  }
  if(VolQ<=0)
  {
     VolQ = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Volume of the domain containing particles should be positive. Set the default value!");
  }
  if( (ParRadius<=0) || (ParRadius>=cbrt(VolQ)) )
  {
     ParRadius = 1.0e-3;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Particle radius should be positive and less than the domain size. Set the default value!");
  }
  if( (ParticleDistance<=0) || (ParticleDistance>=cbrt(VolQ)) )
  {
     ParticleDistance = pow(ParRadius,(2-Kappa)/3);   //ParticleDistance = O(ParticleRadius^(1/3))  
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Distance between neighboring particles should be positive and less than the domain size. Set the default value!");
  }
  if(cnorm<PetscReal,PetscScalar>(OriginalRefractionCoef)<=0)
  {
     OriginalRefractionCoef = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Original refraction should not be zero. Set the default value!");
  }
  if(cnorm<PetscReal,PetscScalar>(DesiredRefractionCoef)<=0)
  {
     DesiredRefractionCoef = -1+I*0.001;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Desired refraction should not be zero. Set the default value!");
  }
  if(Distribution<=0)
  {
     Distribution = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Distribution of particles should be positive. Set the default value!");
  }
  if(TotalParticles<=0)
  {
     TotalParticles = pow(round(1/ParticleDistance),3);   //TotalParticles = O(1/ParticleRadius)  
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Total number of particles should be positive. Set the default value!");
  }
  if( (TotalCubes<=0) || (TotalCubes>TotalParticles) || !is_nth_power(TotalCubes,3) )
  {
     TotalCubes = pow(round(pow(TotalParticles,1.0/3)/10),3);
     if(TotalCubes<=0)
        TotalCubes = 1;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Total number of sub cubes should be positive, cubic and less than total particles. Set the default value!");
  }
  if((TotalCollocationPoints>TotalParticles) || (TotalCollocationPoints<TotalCubes) || !is_nth_power(TotalCollocationPoints,3))
  {
     TotalCollocationPoints = TotalCubes;
     PetscPrintf(PETSC_COMM_WORLD,"\nWARNING: Total number of collocation points should be positive, cubic, and in (total cubes,total particles). Set the default value!");
  }
}

template <class RealVector, class Complex, class Real, class Long, class Integer> 
void Output(Real kappa, Real VolQ, Real ParRadius, Real ParDist, Long TotalParticles, RealVector WaveDirection,\
            Complex OriginalRefractionCoef, Complex DesiredRefractionCoef, Real Distribution, Integer TotalCubes, Complex BoundaryImpedance, Integer& TotalCollocationPoints)
{
    //ostringstream stream;
    string s;

    cout.precision(6);
    cout<< "\n----------------------- SOLVING PARTICLE SCATTERING PROBLEM -----------------------\n";

    cout<< "\nSpeed of wave, v:\t\t\t\t\t";
    cout<< C;
    cout<< "\nFrequency, f:\t\t\t\t\t\t";
    cout<< F;
    cout<< "\nWave number, k:\t\t\t\t\t\t";
    cout<< K;    
    cout<< "\nDirection of plane wave, alpha:\t\t\t\t(";
    cout<< WaveDirection[0]<<", "<<WaveDirection[1]<<", "<<WaveDirection[2]<<")";
    cout<< "\nKappa:\t\t\t\t\t\t\t";
    cout<< kappa;
    cout<< "\nVolume of the domain that contains all particles, |D|:\t";
    cout<< VolQ;
    cout<< "\nOriginal refraction coefficient, n0:\t\t\t";
    s = (cimag(OriginalRefractionCoef)>=0)?"+":"";
    cout<< creal(OriginalRefractionCoef)<<s<<cimag(OriginalRefractionCoef)<<"i";
    cout<< "\nDesired refraction coefficient, n:\t\t\t";
    s = (cimag(DesiredRefractionCoef)>=0)?"+":"";
    cout<< creal(DesiredRefractionCoef)<<s<<cimag(DesiredRefractionCoef)<<"i";    
    cout<< "\nDistribution of particles, N:\t\t\t\t";
    cout<< Distribution;
    cout<< "\nBoundary impedance, zeta:\t\t\t\t";
    s = (cimag(BoundaryImpedance)>=0)?"+":"";
    cout<< creal(BoundaryImpedance)<<s<<cimag(BoundaryImpedance)<<"i"; 
    cout<< "\nRadius of one particle, a:\t\t\t\t";
    cout<< scientific<<ParRadius;
    cout<< "\nDistance between two neighboring particles, d:\t\t";
    cout<< scientific<<ParDist;
    cout<< "\nNumber of particles, M:\t\t\t\t\t";
    cout<< TotalParticles<<" ("<<scientific<<Real(TotalParticles)<<")";
    cout<< "\nNumber of sub cubes after partitioning the domain, P:\t";
    cout<< TotalCubes;
    cout<< "\nTotal collocation points for solving integral equation:\t";
    cout<< TotalCollocationPoints;
    
    cout<<endl;
} 

void checkTime(time_t StartTime, time_t EndTime)
{
    int hour, min, sec, time;

    cout<<"\nStarted on:\t"<<ctime(&StartTime);
    cout<<"Finished on:\t"<<ctime(&EndTime);    	
    time = EndTime-StartTime;
    hour=time/3600;
    time=time%3600;
    min=time/60;
    time=time%60;
    sec=time;
    cout<<"Elapsed time:\t"<<hour<<":"<<min<<":"<<sec;
    cout<<endl;
}

#endif //COMMON_H
