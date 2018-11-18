#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdio.h>

#define h_smooth 1.0e-1//smoothing length
#define k1 1.0 //pressure constant
#define dt 1.0e-4//time step size
#define rho0 1000.0 // reference density
#define m (M_PI*h_smooth*h_smooth*rho0/12.0) //particle mass
#define rigidMassMultiplier 5.0
#define rigidMass m*rigidMassMultiplier
#define nu 2.0e-2 //viscosity coefficient
#define g 9.8//gravitational constant
#define kappa 10.0
#define Ch 15/(14*M_PI*h*h) //normalization constant of cubic spline
#define cs 80.0
//#define dh 0.0001 
#define epsilon 1.0e-8 //small number not to make denominator in gradKernel zero
#define T 50000//time step
#define DAMPTIME 2000//2.5sec for exp2
#define MOTION_START_TIME DAMPTIME+2000//at 3.9sec impact happens for exp2
#define ENLARGEMENT 3.0
#define FLUID_INTERACTION 60.0
#define CONTACT_ANGLE M_PI/4.0
#define HPHILY_INTERACTION FLUID_INTERACTION*(cos(CONTACT_ANGLE))
#define HPHOBY_INTERACTION -FLUID_INTERACTION*(cos(CONTACT_ANGLE)) //negative value
#define IMPACT_VELOCITY 0.0
#define ANGLE_OF_INCIDENT M_PI/2.0
#define BOUNDARY_FORCE 0//if this value is zero calcAccelByBoundaryForce is not called.

#define interval h_smooth
#define MAX_X 30//40
#define MAX_Y 20//20
#define MIN_X -10//-20
#define MIN_Y -5//-5
#define BktLgth 0.3//0.3
#define BktNum (double)(1.0/BktLgth)
#define nBx ((int)((MAX_X-MIN_X)/BktLgth)+2 )
#define nBy ((int)((MAX_Y-MIN_Y)/BktLgth)+2 )
#define nBxy (nBx*nBy)

#endif //_PARAMETERS_H_
