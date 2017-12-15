#ifndef _SPH_H_
#define _SPH_H_

#include <stdio.h>

#define h 0.1//smoothing length
#define k1 1.0 //pressure constant
#define dt 0.5e-4//time step size
#define rho0 1000.0 // reference density
#define m M_PI*h*h*rho0/12.0 //particle mass
#define nu 0.0005 //viscosity coefficient
#define g 9.8//gravitational constant
#define gamm 1.0//surface tension coefficient
#define Ch 15/(14*M_PI*h*h) //normalization constant of cubic spline
#define cs 88.5
//#define dh 0.0001 
#define epsilon 1.0e-5 //small number not to make denominator in gradKernel zero
#define T 170000//time step
#define DAMPTIME 50000//2.5sec for exp2
#define ROTSTARTTIME DAMPTIME+1//at 3.9sec impact happens for exp2
#define FLUID_INTERACTION 1.0
#define HPHILY_INTERACTION 0.5
#define HPHOBY_INTERACTION -0.5//negative value
#define IMPACT_VELOCITY 5.0
#define ANGLE_OF_INCIDENT M_PI/2.0

#define interval 0.1
#define MAX_X 60
#define MAX_Y 60
#define MIN_X -5  
#define MIN_Y -5
#define BktLgth 0.3
#define BktNum 1.0/BktLgth
#define nBx ((int)((MAX_X-MIN_X)/BktLgth)+2 )
#define nBy ((int)((MAX_Y-MIN_Y)/BktLgth)+2 )
#define nBxy (nBx*nBy)


typedef struct{
  int inRegion;
  short color;
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
  double rho; //density
  double p; //pressure
  double mass;
}Particle_State;

double cubicSpline1(double q);
double cubicSpline2(double q);
double kernel(Particle_State p1, Particle_State p2);
double gradKernel(Particle_State p1, Particle_State p2, int x_or_y);
double Laplacian(Particle_State p1, Particle_State p2);
void calcDensity(Particle_State p[], int bfst[], int nxt[]);
void calcPressure(Particle_State p[]);
void initializeAccel(Particle_State p[]);
void calcAccelByExternalForces(Particle_State p[]);
void calcAccelByPressure(Particle_State p[], int bfst[], int nxt[]);
void calcAccelByViscosity(Particle_State p[], int bfst[], int nxt[], int time);
double surfaceTensionCoefficient(double r);
void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int nxt[]);
void calcInterfacialForce(Particle_State p[], int bfst[], int nxt[]);
double boundaryGamma(Particle_State p1, Particle_State p2);
void calcAccelByBoundaryForce(Particle_State p[], int bfst[], int nxt[]);
double adhesionCoefficient(Particle_State p1, Particle_State p2);
void calcAccelByAdhesion(Particle_State p[], int bfst[], int nxt[]);
void rotateRigidBody(Particle_State p[], double angVel);
void rigidBodyCorrection(Particle_State p[], FILE *fp, int time, double com[]);
void leapfrogStart(Particle_State p[]);
void leapfrogStep(Particle_State p[], int time);

void initialization(Particle_State p[], int particleNumber);
int fluidParticles(Particle_State p[]);
int wallParticles(Particle_State p[]);
int obstacleBoundaryParticles(Particle_State p[]);
void setInitialVelocity(Particle_State p[]);

void allocateBucket(int **bfst, int **blst, int **nxt);
void checkParticle(Particle_State p[]);
void makeBucket(int *bfst, int *blst, int*nxt, Particle_State p[]);
void freeBucket(int *bfst, int *blst, int*nxt);
double gradSpikey(Particle_State p1, Particle_State p2, int axis);
double poly6(Particle_State p1, Particle_State p2);

void getSourceImageName(FILE *fp, char srcName[]);
double calcRadius(Particle_State p[]);
void makeFileNamePrefix(char fNamePrefix[], char srcName[], double angVel, double spinParam);
void printParticles(Particle_State p[], FILE *fp);
void percentage(int time, int *countPer);
void printFluidParticles(Particle_State p[], FILE *fp);
void printFluidPositions(Particle_State p[], FILE *fp);
void printBoundaryParticles(Particle_State p[], FILE *fp);
void printBoundaryPositions(Particle_State p[], FILE *fp);
void printObstacleParticles(Particle_State p[], FILE *fp);
void printObstaclePositions(Particle_State p[], FILE *fp);
void printParticlesAroundObstacle(Particle_State p[], FILE *fp, double com[]);
void printParameters(FILE *fp, double angVel, char srcName[], char date[], double spinParam);
void tipPosition(Particle_State p[], int time, FILE *tip);
void getCalculationRegion(double range[], Particle_State p[]);
void makePltFile(char *srcName, Particle_State p[], char *fileNamePrefix);
#endif //_SPH
