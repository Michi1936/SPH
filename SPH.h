#ifndef _SPH_H_
#define _SPH_H_

#include <stdio.h>

#define h 1.0e-1//smoothing length
#define k1 1.0 //pressure constant
#define dt 0.5e-4//time step size
#define rho0 1000.0 // reference density
#define m M_PI*h*h*rho0/12.0 //particle mass
#define rigidMass m*2.0
#define nu 1.0e-3 //viscosity coefficient
#define g 9.8//gravitational constant
#define kappa 0.0
#define Ch 15/(14*M_PI*h*h) //normalization constant of cubic spline
#define cs 100
//#define dh 0.0001 
#define epsilon 1.0e-8 //small number not to make denominator in gradKernel zero
#define T 25000//time step
#define DAMPTIME 4000//2.5sec for exp2
#define MOTION_START_TIME DAMPTIME+1000//at 3.9sec impact happens for exp2
#define FLUID_INTERACTION 0.0
#define HPHILY_INTERACTION FLUID_INTERACTION/2.0
#define HPHOBY_INTERACTION -FLUID_INTERACTION/2.0//negative value
#define IMPACT_VELOCITY 10.0
#define ANGLE_OF_INCIDENT M_PI/2.0
#define BOUNDARY_FORCE 1//if this value is zero calcAccelByBoundaryForce is not called.

#define interval 0.1
#define MAX_X 50
#define MAX_Y 50
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

typedef struct{
  double prepx;
  double prepy;
}RigidPreValue;

double cubicSpline1(double q);
double cubicSpline2(double q);
double kernel(Particle_State p1, Particle_State p2);
double gradKernel(Particle_State p1, Particle_State p2, int x_or_y);
double Laplacian(Particle_State p1, Particle_State p2);
void calcPsi(Particle_State p[], double Psi[], int bfst[], int nxt[]);
void AkinciCalcDensity(Particle_State p[], double Psi[], int bfst[], int nxt[]);
void calcPressure(Particle_State p[]);
void initializeAccel(Particle_State p[]);
void calcAccelByExternalForces(Particle_State p[]);
void AkinciCalcAccelByPressure(Particle_State p[], double Psi[], int bfst[], int nxt[]);
void AkinciCalcAccelByViscosity(Particle_State p[], double Psi[], int bfst[], int nxt[], int time);
double surfaceTensionCoefficient(double r);
void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int nxt[]);
void calcInterfacialForce(Particle_State p[], int bfst[], int nxt[]);
double boundaryGamma(Particle_State p1, Particle_State p2);
void calcAccelByBoundaryForce(Particle_State p[], int bfst[], int nxt[]);
double adhesionCoefficient(Particle_State p1, Particle_State p2);
void calcAccelByAdhesion(Particle_State p[], int bfst[], int nxt[]);
void rotateRigidBody(Particle_State p[], RigidPreValue rig[], double angVel);
//void rigidBodyCorrection(Particle_State p[], RigidPreValue rig[],FILE *fp, int time, double com[]);
void leapfrogStart(Particle_State p[], RigidPreValue rig[]);
void leapfrogStep(Particle_State p[], RigidPreValue rig[], int time);
void rigidBodyTimeIntegration(Particle_State p[], double *omega, FILE *fp, int time);

void initialization(Particle_State p[], RigidPreValue rig[]);
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
void getMaxVelocity(Particle_State p[], FILE *fp, int time);
void getCalculationRegion(double range[], Particle_State p[]);
void makePltFile(char *srcName, Particle_State p[], char *fileNamePrefix);

#endif //_SPH
