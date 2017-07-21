#ifndef _SPH_H_
#define _SPH_H_

#define m 0.003 //particle mass
#define h 0.07 //smoothing length
#define k1 15.0 //pressure constant
#define dt 6.0e-3 //time step size
#define rho0 3.0 // reference density
#define nu 0.1 //viscosity coefficient
#define g 9.8//gravitational constant
#define gamm 0//surface tension coefficient
#define Ch 15/(14*M_PI*h*h) //normalization constant of cubic spline
#define dh 0.0001 
#define epsilon 1.0e-5 //small number not to make denominator in gradKernel zero
#define T 400//time step

#define MAX_X 5   
#define MAX_Y 4  
#define MIN_X -5   
#define MIN_Y -10
#define BktLgth 0.3
#define BktNum 1.0/BktLgth
#define nBx ((int)((MAX_X-MIN_X)/BktLgth)+2 )
#define nBy ((int)((MAX_Y-MIN_Y)/BktLgth)+2 )
#define nBxy (nBx*nBy)


typedef struct{
  int inRegion;
  double px; 
  double py;
  double vx;
  double vy;
  double vxh;
  double vyh;
  double ax;
  double ay;
  double rho; //density
  double p; //pressure
}Particle_State;

double cubicSpline1(double q);
double cubicSpline2(double q);
double cohesion(double r);
double kernel(Particle_State p1, Particle_State p2);
double gradKernel(Particle_State p1, Particle_State p2, int x_or_y);
double Laplacian(Particle_State p1, Particle_State p2);
void calcDensity(Particle_State p[], int bfst[], int blst[], int nxt[]);
void calcPressure(Particle_State p[]);
void calcAcceleration(Particle_State p[], int bfst[], int blst[], int nxt[]);
void timeDevelopment(Particle_State p[]);
void leapfrogStart(Particle_State p[]);
void leapfrogStep(Particle_State p[]);
void initialization(Particle_State p[], int particleNumber);
int fluidParticles(Particle_State p[]);
int wallParticles(Particle_State p[]);
void obstacleBoundaryParticles(Particle_State obp[]);
void allocateBucket(int **bfst, int **blst, int **nxt);
void checkParticle(Particle_State p[]);
void makeBucket(int *bfst, int *blst, int*nxt, Particle_State p[]);
void freeBucket(int *bfst, int *blst, int*nxt);
void normalizePressure(Particle_State p[]);

#endif //_SPH
