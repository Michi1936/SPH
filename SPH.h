#ifndef _SPH_H_
#define _SPH_H_

#define k1 1.0 //pressure constant
#define rho0 1000.0 // reference density
#define h 0.1//smoothing length
#define k1 1.0 //pressure constant
#define dt 2.0e-4//time step size
#define rho0 1000.0 // reference density
#define m M_PI*h*h*rho0/12.0 //particle mass
#define nu 0.001 //viscosity coefficient
#define g 9.8//gravitational constant
#define gamm 1.0//surface tension coefficient
#define Ch 15/(14*M_PI*h*h) //normalization constant of cubic spline
#define cs 88.5
//#define dh 0.0001 
#define epsilon 1.0e-5 //small number not to make denominator in gradKernel zero
#define T 10000//time step



#define interval 0.1

#define MAX_X 25   
#define MAX_Y 25
#define MIN_X -10  
#define MIN_Y -10
#define BktNum 1.0/BktLgth
#define BktLgth 0.3
#define nBx ((int)((MAX_X-MIN_X)/BktLgth)+2 )
#define nBy ((int)((MAX_Y-MIN_Y)/BktLgth)+2 )
#define nBxy (nBx*nBy)


typedef struct{
  int inRegion;
  double px; 
  double py;
  double prepx;
  double prepy;
  double vx;
  double vy;
  double vxh;
  double vyh;
  double ax;
  double ay;
  double mass;//mass of particle
  double rho; //density
  double p; //pressure
  double Theta;//inverse of the particle volume
  double mu;
}Particle_State;

double cubicSpline1(double q);
double cubicSpline2(double q);
double kernel(Particle_State p1, Particle_State p2);
double gradKernel(Particle_State p1, Particle_State p2, int x_or_y);
double Laplacian(Particle_State p1, Particle_State p2);
void calcInverseParticleVolume(Particle_State p[], int bfst[], int nxt[]);
void calcDensity(Particle_State p[], int bfst[], int blst[], int nxt[]);
void calcPressure(Particle_State p[]);
void initializeAccel(Particle_State p[]);
void calcAccelByExternalForces(Particle_State p[], int bfst[], int blst[], int nxt[]);
void calcAccelByPressure(Particle_State p[], int bfst[], int nxt[]);
void calcAccelByViscosity(Particle_State p[], int bfst[], int blst[], int nxt[]);
double surfaceTensionCoefficient(double r);
void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int blst[], int nxt[]);
double boundaryGamma(Particle_State p1, Particle_State p2);
void calcAccelByBoundaryForce(Particle_State p[], int bfst[], int nxt[]);
double adhesionCoefficient(Particle_State p1, Particle_State p2);
void calcAccelByAdhesion(Particle_State p[], int bfst[], int nxt[]);
void rigidBodyCorrection(Particle_State p[], int bfst[], int nxt[]);
void leapfrogStart(Particle_State p[]);
void leapfrogStep(Particle_State p[]);

void initialization(Particle_State p[], int particleNumber);
int fluidParticles(Particle_State p[]);
int airParticles(Particle_State p[]);
int wallParticles(Particle_State p[]);
int obstacleBoundaryParticles(Particle_State p[]);

void allocateBucket(int **bfst, int **blst, int **nxt);
void checkParticle(Particle_State p[]);
void makeBucket(int *bfst, int *blst, int*nxt, Particle_State p[]);
void freeBucket(int *bfst, int *blst, int*nxt);


#endif //_SPH
