#ifndef _SPH_H_
#define _SPH_H_

#include <stdio.h>

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
  double angle;
  double omega;
}RigidBodyValues;

double cubicSpline1(double q);
double cubicSpline2(double q);
double kernel(Particle_State p1, Particle_State p2);
double gradKernel(Particle_State p1, Particle_State p2, int x_or_y);
double Laplacian(Particle_State p1, Particle_State p2);
void calcDensity(Particle_State p[], int bfst[], int nxt[]);
void calcPressure(Particle_State p[]);
void initializeAccel(Particle_State p[]);
void calcAccelByExternalForces(Particle_State p[]);
double surfaceTensionCoefficient(double r);
void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int nxt[]);
void calcInterfacialForce(Particle_State p[], int bfst[], int nxt[]);
double boundaryGamma(Particle_State p1, Particle_State p2);
void calcAccelByBoundaryForce(Particle_State p[], int bfst[], int nxt[]);
void setInitialVelocity(Particle_State p[]);
void rotateRigidBody(Particle_State p[], RigidBodyValues rig, double angVel);

void set_boundaryType(int boundaryType[]);
void calcPsi(Particle_State p[], double Psi[], int bfst[], int nxt[], int rigidNum[]);
void AkinciCalcDensity(Particle_State p[], double Psi[], int bfst[], int nxt[]);
void AkinciCalcAccelByPressure(Particle_State p[], double Psi[], int bfst[], int nxt[]);
void AkinciCalcAccelByViscosity(Particle_State p[], double Psi[], int bfst[], int nxt[], int time);
void EulerCromerTimeIntegration(Particle_State p[]);
void rigidBodyTimeIntegration(Particle_State p[], RigidBodyValues rigV, FILE *fp, int time);

void initialization(Particle_State p[], RigidBodyValues rig);
int fluidParticles(Particle_State p[]);
int wallParticles(Particle_State p[]);
int obstacleBoundaryParticles(Particle_State p[]);


void allocateBucket(int **bfst, int **blst, int **nxt);
void checkParticle(Particle_State p[]);
void makeBucket(int *bfst, int *blst, int*nxt, Particle_State p[]);
void freeBucket(int *bfst, int *blst, int*nxt);
double gradSpikey(Particle_State p1, Particle_State p2, int axis);
double poly6(Particle_State p1, Particle_State p2);

void getSourceImageName(FILE *fp, char srcName[]);
double calcRadius(Particle_State p[]);
void makeFileNamePrefix(char fNamePrefix[], char srcName[], double angVel, double spinParam);
void openDatFile(FILE **fp, char type[], char srcName[], char prefix[]);
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

void getMaxVelocity(Particle_State p[], FILE *fp, int time);
void getCalculationRegion(double range[], Particle_State p[]);
void makePltFile(char *srcName, Particle_State p[], char *fileNamePrefix);

#endif //_SPH
