#include "common.h"



///////////////////////////////////////////////////////////////////
template < class RealVector, class ComplexVector, class Complex, class Real, class Long, class Integer >
class Scattering3DP 
{
    public:
        Integer        NumParticlePerSide;
        Integer        NumParticlePerPlane;
        Integer        NumCubePerSide;
        Integer        NumCubePerPlane;
        Integer        NumParticlePerSmallCubeSide;
        Real           DomainSize;
        Real           SmallCubeSize;
        Real           ParticleRadius;
        Real           ParticleDistance;
        Real           Kappa;
        Long           TotalParticles;
        Integer        TotalCubes;
        Integer        TotalCollocationPoints;
        RealVector     WaveDirection;
        //ComplexVector  u0; //Initial field
        Complex        OriginalRefractionCoef;
        Complex        DesiredRefractionCoef;
        Real           Distribution;
        Real           VolQ;  //Volume of the cube that contains all particles
        Real           SmallCubeVol;
        Complex        p;     
        Complex        h;
        Complex        BoundaryImpedance;

    private:
        Complex        ha2k;  //coeficient of entries in Matrix A
        Complex        ha2kPI4;  //coeficient of entries in Matrix A
        Real           iNumCubePerPlane; //inverse(NumCubePerPlane)
        Real           iNumCubePerSide;  //inverse(NumCubePerSide)
        RealVector     x;     //Particles positions
        RealVector     y;
        RealVector     z;

    public:
        Scattering3DP()
        {

        }

        ///////////////////////////////////////////////////////////////////

        virtual ~Scattering3DP()
        {

        }

        ///////////////////////////////////////////////////////////////////
        void Input(Real ParticleRadius, Real Kappa, RealVector WaveDirection, Real ParticleDistance, Long TotalParticles)
        {
            this->ParticleRadius = ParticleRadius;
            this->Kappa = Kappa;
            this->WaveDirection = WaveDirection;
            this->ParticleDistance = ParticleDistance;
            this->TotalParticles = TotalParticles;
        }

        ///////////////////////////////////////////////////////////////////
        void Input(Real ParticleRadius, Real Kappa, RealVector WaveDirection, Real ParticleDistance, Long TotalParticles,\
                   Complex OriginalRefractionCoef, Complex DesiredRefractionCoef, Real Distribution, Real VolQ, Integer TotalCubes, \
                   Integer TotalCollocationPoints)
        {
            this->ParticleRadius = ParticleRadius;
            this->Kappa = Kappa;
            this->WaveDirection = WaveDirection;
            this->ParticleDistance = ParticleDistance;
            this->TotalParticles = TotalParticles;
            this->OriginalRefractionCoef = OriginalRefractionCoef;
            this->DesiredRefractionCoef = DesiredRefractionCoef;
            this->Distribution = Distribution;
            this->VolQ = VolQ; 
            this->TotalCubes = TotalCubes;  
            this->TotalCollocationPoints = TotalCollocationPoints;      
        }


        ///////////////////////////////////////////////////////////////////

        void Init()
        {            
            // Number of particles on a side of a cube of size 1
            NumParticlePerSide = round(pow(TotalParticles,1.0/3));
            NumParticlePerPlane = pow(NumParticlePerSide,2);
            NumCubePerSide = ceil(pow(TotalCubes,1.0/3));
            NumCubePerPlane = pow(NumCubePerSide,2);
            DomainSize = (NumParticlePerSide-1)*ParticleDistance;
            SmallCubeSize = DomainSize/NumCubePerSide;
            SmallCubeVol = pow(SmallCubeSize,3);
            NumParticlePerSmallCubeSide = floor(NumParticlePerSide/NumCubePerSide);
 
            //Performance tune:
            iNumCubePerPlane = 1.0/NumCubePerPlane;
            iNumCubePerSide = 1.0/NumCubePerSide;

            // Initial field satisfies Helmholtz equation in R^3
            //u0 = InitField();

            //Cubes positions            
            UniformDistributeCubes();
                      
            p = pow(K,2)*(cpow(OriginalRefractionCoef,2) - cpow(DesiredRefractionCoef,2));
            Real h1 = creal(p)/(PI4*Distribution);
            Real h2 = cimag(p)/(PI4*Distribution);
            h = (h1+I*h2);
            ha2k = h*Distribution*SmallCubeVol;
            ha2kPI4 = ha2k*PI4;
            BoundaryImpedance = h/pow(ParticleRadius,Kappa);
        }

        ///////////////////////////////////////////////////////////////////

        void UniformDistributeParticles()
        {
            // Set the position for each particle (uniformly distributed)
            Real x0,y0,z0,t;

            //Set the cube centered at the origin
            t = pow(VolQ,1.0/3)/2;
            x0 = -t;
            y0 = -t;
            z0 = -t;

            // The first particle [x1,y1,z1] is at the left bottom corner of the
            // cube and is called particle number 1.
            x = RealVector(NumParticlePerSide);
            y = RealVector(NumParticlePerSide);
            z = RealVector(NumParticlePerSide);
            
            for(Integer s = 0; s < NumParticlePerSide; s++)
            {
                t = ParticleDistance*s;
                x[s] = x0 + t;
                y[s] = y0 + t;
                z[s] = z0 + t;
            }
        }

        ///////////////////////////////////////////////////////////////////

        void UniformDistributeCubes()
        {
            // Set the position for each cube (uniformly distributed)
            Real x0,y0,z0,t;

            //Set the cube centered at the origin
            t = pow(VolQ,1.0/3)/2;
            x0 = -t+SmallCubeSize/2;
            y0 = x0;
            z0 = x0;

            // The first small cube [x1,y1,z1] is at the left bottom corner of the
            // cube and is called cube number 1.
            x = RealVector(NumCubePerSide);
            y = RealVector(NumCubePerSide);
            z = RealVector(NumCubePerSide);
            
            for(Integer s = 0; s < NumCubePerSide; s++)
            {
                t = SmallCubeSize*s;
                x[s] = x0 + t;
                y[s] = y0 + t;
                z[s] = z0 + t;
            }
        }

        ///////////////////////////////////////////////////////////////////
/*
        inline const Real DistributionFunc(const Long& CubeNumber)
        {
            return 1;
        }

        ///////////////////////////////////////////////////////////////////

        inline const Complex DesiredRefractionCoefFunc(const Long& CubeNumber)
        {
            return sqrt(0.2);
        }

        ///////////////////////////////////////////////////////////////////

        inline const Complex OriginalRefractionCoefFunc(const Long& CubeNumber)
        {
            return 1;
        }

        ///////////////////////////////////////////////////////////////////

        inline const Complex p_Func(const Long& CubeNumber)
        {
            return cpow(K,2)*(cpow(OriginalRefractionCoefFunc(CubeNumber),2) - cpow(DesiredRefractionCoefFunc(CubeNumber),2));
        }

        ///////////////////////////////////////////////////////////////////

        inline const Complex h_Func(const Long& CubeNumber)
        {
            Complex p = p_Func(CubeNumber);
            Real h1 = creal(p)/(PI4*DistributionFunc(CubeNumber));
            Real h2 = cimag(p)/(PI4*DistributionFunc(CubeNumber));

            return (h1+I*h2);
        }
*/
        ///////////////////////////////////////////////////////////////////
        //Uniform Distribution
        inline const RealVector Particle2Position(Long s)
        {
            // Return the position in the 3D cube of particle s

            // The first particle [x1,y1,z1] is at the left, bottom corner of the
            // cube and is called particle number 1. The next one will be on the same
            // row, go to the right. When finishing the first line, go to the second line
            // and start at the first column again. When finishing the first plane, move
            // up.

            // [x1,x2,x3] is an array index
            Integer x3 = floor(s/(NumParticlePerPlane)); // Find the plane where the particle s is on
            Integer x2 = s%NumParticlePerSide;
            Integer t = s%(NumParticlePerPlane);
            Integer x1 = floor(t/NumParticlePerSide);

            RealVector v(3);
            v[0] = x[x1];
            v[1] = y[x2];
            v[2] = z[x3];

            return v;
        }

        ///////////////////////////////////////////////////////////////////
        //Uniform Distribution - this function increases performance
        inline const void Particle2Position(Long s,Real& xs,Real& ys,Real& zs)
        {
            // Return the position in the 3D cube of particle s

            // The first particle [x1,y1,z1] is at the left, bottom corner of the
            // cube and is called particle number 1. The next one will be on the same
            // row, go to the right. When finishing the first line, go to the second line
            // and start at the first column again. When finishing the first plane, move
            // up.

            // [x1,x2,x3] is an array index
            Integer x1,x2,x3,t;
            x3 = floor(s/NumParticlePerPlane); // Find the plane where the particle s is on
            x2 = s%NumParticlePerSide;
            t = s%(NumParticlePerPlane);
            x1 = floor(t/NumParticlePerSide);
                     
            xs = x[x1];
            ys = y[x2];
            zs = z[x3];
        }

        ///////////////////////////////////////////////////////////////////

        //Uniform Distribution
        inline const RealVector Cube2Position(Integer s)
        {
            // Return the position in 3D of cube s

            // The first small cube [x1,y1,z1] is at the left, bottom corner of the big
            // cube and is called small cube number 1. The next one will be on the same
            // row, go to the right. When finishing the first line, go to the second line
            // and start at the first column again. When finishing the first plane, move
            // up.

            // [x1,x2,x3] is an array index
            Integer x1,x2,x3,t;
            x3 = floor(s/NumCubePerPlane); // Find the plane where the cube s is on
            x2 = s%NumCubePerSide;
            t = s%(NumCubePerPlane);
            x1 = floor(t/NumCubePerSide);

            RealVector v(3);
            v[0] = x[x1];
            v[1] = y[x2];
            v[2] = z[x3];

            return v;
        }

        ///////////////////////////////////////////////////////////////////


        //Uniform Distribution - this function increases performance
        inline const void Cube2Position(Integer s,Real& xs,Real& ys,Real& zs)
        {
            // Return the position in 3D of cube s

            // The first small cube [x1,y1,z1] is at the left, bottom corner of the big
            // cube and is called small cube number 1. The next one will be on the same
            // row, go to the right. When finishing the first line, go to the second line
            // and start at the first column again. When finishing the first plane, move
            // up.

            // [x1,x2,x3] is an array index
            Integer x1,x2,x3,t;
            x3 = floor(s*iNumCubePerPlane); // Find the plane where the cube s is on
            x2 = s%NumCubePerSide;
            t = s%(NumCubePerPlane);
            x1 = floor(t*iNumCubePerSide);
                     
            xs = x[x1];
            ys = y[x2];
            zs = z[x3];
        }

        ///////////////////////////////////////////////////////////////////

        Integer FindCube(Integer ParticleNumber)
        {
            //Find cube number that contains particle ParticleNumber        
        
            // [x1,x2,x3] is an array index of particle ParticleNumber
            Integer x1,x2,x3,t,CubeNumber;

            x3 = floor(ParticleNumber/NumParticlePerPlane); // Find the plane where the particle ParticleNumber is on
            x2 = ParticleNumber%NumParticlePerSide;
            t = ParticleNumber%(NumParticlePerPlane);
            x1 = floor(t/NumParticlePerSide);

            //Find the index of the small cube containing particle ParticleNumber            
            x1 = floor(x1/NumParticlePerSmallCubeSide); 
            x2 = floor(x2/NumParticlePerSmallCubeSide); 
            x3 = floor(x3/NumParticlePerSmallCubeSide); 

            if(x1>=NumCubePerSide)
                 x1 = NumCubePerSide - 1;
            if(x2>=NumCubePerSide)
                 x2 = NumCubePerSide - 1;
            if(x3>=NumCubePerSide)
                 x3 = NumCubePerSide - 1;

            CubeNumber = x2 + x1*NumCubePerSide + x3*NumCubePerPlane;

            return CubeNumber;                 
        }

        ///////////////////////////////////////////////////////////////////

        Integer FindCubeOfCollocation(Integer CollocationPoint)
        {
            //Find cube number that contains particle CollocationPoint
        
            // [x1,x2,x3] is an array index of the CollocationPoint
            Integer x1,x2,x3,t,CubeNumber;

            Integer NumCollocationPointsPerSide = round(pow(TotalCollocationPoints,1.0/3));
            Integer NumCollocationPointsPerPlane = pow(NumCollocationPointsPerSide,2);
            Integer NumCollocationPointsPerSmallCubeSide = floor(NumCollocationPointsPerSide/NumCubePerSide);

            x3 = floor(CollocationPoint/NumCollocationPointsPerPlane); // Find the plane where the CollocationPoint is on
            x2 = CollocationPoint%NumCollocationPointsPerSide;
            t = CollocationPoint%(NumCollocationPointsPerPlane);
            x1 = floor(t/NumCollocationPointsPerSide);

            //Find the index of the small cube containing the CollocationPoint            
            x1 = floor(x1/NumCollocationPointsPerSmallCubeSide); 
            x2 = floor(x2/NumCollocationPointsPerSmallCubeSide); 
            x3 = floor(x3/NumCollocationPointsPerSmallCubeSide); 

            if(x1>=NumCubePerSide)
                 x1 = NumCubePerSide - 1;
            if(x2>=NumCubePerSide)
                 x2 = NumCubePerSide - 1;
            if(x3>=NumCubePerSide)
                 x3 = NumCubePerSide - 1;

            CubeNumber = x2 + x1*NumCubePerSide + x3*NumCubePerPlane;

            return CubeNumber;                 
        }


        ///////////////////////////////////////////////////////////////////


        inline const Complex Green(Long s, Long t)
        {
            // Create a Green function in R^3, s!=t            

            // Distance from cube s to cube t in R^3 
            /*
            RealVector spos = Cube2Position(s);
            RealVector tpos = Cube2Position(t);
            RealVector d(3);
            d[0] = spos[0]-tpos[0];
            d[1] = spos[1]-tpos[1];
            d[2] = spos[2]-tpos[2];
            Real st = Norm<RealVector,Real>(d, 2);
            */   
            //For performance:
            Real xs,ys,zs,xt,yt,zt;
            Cube2Position(s,xs,ys,zs);
            Cube2Position(t,xt,yt,zt);          
            Real st = sqrt(pow(xs-xt,2)+pow(ys-yt,2)+pow(zs-zt,2)); 

            if(st<=ZERO)
                return 0;

            Complex G = cexp(I*K*st)/(PI4*st);
	     //cout<<"\nnorm(G["<<s<<","<<t<<"])="<<cnorm<Real, Complex>(G);

            return G;
        }

        ///////////////////////////////////////////////////////////////////

        inline ComplexVector InitField()
        {
            // Create an inittial field u0 satisfying Helmholtz equation in R^3

            ComplexVector u0(TotalParticles);
            RealVector parpos;

            for(Long s = 0; s < TotalParticles; s++)
            {
                parpos = Cube2Position(s);
                u0[s] = cexp(I*K*Dot<Real,RealVector>(WaveDirection,parpos));
            }

            return u0;
        }

        ///////////////////////////////////////////////////////////////////

        inline const Complex InitField(Long s)
        {
            // Create an inittial field u0 satisfying Helmholtz equation in R^3
            /*
            RealVector parpos = Cube2Position(s);
            Complex u0 = cexp(I*K*Dot<Real,RealVector>(WaveDirection,parpos));            
            */
            //For performance:
            Real xs,ys,zs;
            Cube2Position(s,xs,ys,zs);
            Complex u0 = cexp(I*K*(WaveDirection[0]*xs + WaveDirection[1]*ys + WaveDirection[2]*zs)); 

            return u0;
        }

        ///////////////////////////////////////////////////////////////////
/*
        inline const ComplexVector AuFunc(const ComplexVector u)
        {
            // Compute A*u

            ComplexVector Au(TotalCubes);

            for(Long s = 0; s < TotalCubes; s++)
            {
                for(Long t = 0; t < TotalCubes; t++)
                {
                    if (s!=t)
                        Au(s) = Au(s) + (Green(s,t)*h_Func(t)*DistributionFunc(t)*SmallCubeVol*PI4)*u(t);
                    else
                        Au(s) = Au(s) + u(t);

                    //Au(s) += CoefMat(s,t)*u(t);
                }
            }

            return Au;
        }
*/

        ///////////////////////////////////////////////////////////////////
        //Matrix A in Ax=u0
        inline const Complex CoefMat(Long i, Long j)
        {
            // Generate value for entry A(i,j) in Au=u0

            if (i==j)
               return 1;

            return (Green(i,j)*ha2kPI4);
        }

        ///////////////////////////////////////////////////////////////////
        //Matrix A in Ax=u0, Performance tune
        inline const Complex CoefMatFast(Long s, Long t)
        {
            // Generate value for entry A(i,j) in Au=u0

            if (s==t)
               return 1;

            Real xs,ys,zs,xt,yt,zt;
            Cube2Position(s,xs,ys,zs);
            Cube2Position(t,xt,yt,zt);          
            Real st = sqrt(pow(xs-xt,2)+pow(ys-yt,2)+pow(zs-zt,2)); 

            if(st<=ZERO)
                return 0;

            //Green function*PI4:
            Complex G = cexp(I*K*st)/st;

            return (G*ha2k);
        }
};

