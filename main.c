#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"Parameters.h"
#include"numbers.h"


int main(int argc, char *argv[])
{
  Particle_State a[N];
  RigidPreValue rig[OBP];
  double com[2];
  int *bfst, *blst, *nxt;
  int countPer=0;
  double spinParam=0;
  int i;
  FILE *data;
  FILE *plot;
  FILE *partPlot;
  FILE *parameters;
  FILE *numbers;
  FILE *rigidBody;
  FILE *velocity;
  char srcName[256];
  char fileNamePrefix[256];
  char date[256];
  char type[32];
  double angVel=0;

  time_t t=time(NULL);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("\n\nCalculation started:%s\n", date);

  //defining anglar velocity
  if(argc==2){
    angVel=atof(argv[1]);
  }

  //placing particles
  initialization(a, rig);
  fluidParticles(a);
  wallParticles(a);
  obstacleBoundaryParticles(a);
  fprintf(stderr, "Particles are placed.\n");

  //allocating bucket
  allocateBucket(&bfst, &blst, &nxt);
  fprintf(stderr,"%f %f %d %d %d Bucket allocated\n", BktLgth, BktNum, nBx, nBy, nBxy);
  checkParticle(a);  
  makeBucket(bfst, blst, nxt, a);
  fprintf(stderr,"Bucket is made.\n\n");

  spinParam=fabs(calcRadius(a)*angVel/(IMPACT_VELOCITY+epsilon));
  fprintf(stderr, "spin parameter:%f\n", spinParam);
  
  //open numbers.h
  numbers=fopen("numbers.h","r");
  if(numbers==NULL){
    printf("numbers.h cannot be opened!\n");
    exit(EXIT_FAILURE);
  }

  getSourceImageName(numbers, srcName);
  printf("%s.png\n",srcName);
  makeFileNamePrefix(fileNamePrefix, srcName, angVel, spinParam);
  fprintf(stderr, "\n\nprefix and suffix test\n./Source_%s/%s_plot.dat\n", srcName, fileNamePrefix);

  //opening .dat files
  sprintf(type, "plot");
  openDatFile(&plot, type, srcName, fileNamePrefix);
  sprintf(type, "data");
  openDatFile(&data, type, srcName, fileNamePrefix);
  sprintf(type, "parameters");
  openDatFile(&parameters, type, srcName, fileNamePrefix);
  sprintf(type, "rigidBody");
  openDatFile(&rigidBody, type, srcName, fileNamePrefix);
  sprintf(type, "partPlot");
  openDatFile(&partPlot, type, srcName, fileNamePrefix);
  sprintf(type, "maxVelocity");
  openDatFile(&velocity, type, srcName, fileNamePrefix);

  //printint parameters-------------------
  printParameters(parameters, angVel, srcName, date, spinParam);

  //calculating initial state
  calcDensity(a, bfst, nxt);
  calcPressure(a);
  initializeAccel(a);
  calcAccelByExternalForces(a);
  calcAccelByPressure(a,bfst, nxt);
  calcAccelByViscosity(a,bfst, nxt,0);
  //  calcAccelBySurfaceTension(a, bfst, nxt);

  if(FLUID_INTERACTION>epsilon){
    calcInterfacialForce(a, bfst, nxt);
  }
  if(BOUNDARY_FORCE==1){
    calcAccelByBoundaryForce(a, bfst, nxt);
  }
  
  //printParticles(a,data);//Here shows parameters at t=0
  printBoundaryParticles(a, data);
  printFluidParticles(a, data);
  printObstacleParticles(a, data);

  printBoundaryPositions(a, plot);
  printFluidPositions(a, plot);
  printObstaclePositions(a,plot);

  for(i=FLP+BP; i<N; i++){
    com[0]+=a[i].px/OBP;
    com[1]+=a[i].py/OBP;
  }

  printBoundaryPositions(a, partPlot);
  printParticlesAroundObstacle(a, partPlot, com);

  if(MOTION_START_TIME==0){
    fprintf(stderr, "At 0 rigid body is rotated\n");
    rotateRigidBody(a, rig, angVel);
  }
  
  //time development
  for(i=1; i<=T; i++){
    if(i==1){
      leapfrogStart(a, rig);
    }else{
      leapfrogStep(a,rig,i);
      if(i==MOTION_START_TIME){
        fprintf(stderr, "At %d rigid body is rotated\n", i);
        rotateRigidBody(a, rig, angVel);
        setInitialVelocity(a);
        fprintf(stderr, "Rigid Body Velocity is set\n");
      }
    }

    rigidBodyCorrection(a, rig, rigidBody, i, com);
    velocityCorrection(a, bfst, nxt);   
    checkParticle(a);
    makeBucket(bfst, blst, nxt, a);

    calcDensity(a, bfst, nxt);
    calcPressure(a);
    initializeAccel(a);
    calcAccelByExternalForces(a);
    calcAccelByPressure(a,bfst, nxt);
    calcAccelByViscosity(a,bfst, nxt,i);
<<<<<<< HEAD
    //    calcAccelBySurfaceTension(a, bfst, nxt);

    calcAccelBySurfaceTension(a, bfst, nxt);

>>>>>>> master
    if(FLUID_INTERACTION>epsilon){
      calcInterfacialForce(a, bfst, nxt);
    }
    if(BOUNDARY_FORCE==1){
      calcAccelByBoundaryForce(a, bfst, nxt);
    }
    if(i%100==0){
      printFluidParticles(a, data);//here show paremeters at t=(i*dt)
      printObstacleParticles(a, data);
      printFluidPositions(a,plot);
      printObstaclePositions(a,plot);
      printParticlesAroundObstacle(a, partPlot, com);
    }
    percentage(i, &countPer);
    getMaxVelocity(a, velocity, i);
  }

  fprintf(stderr,"\n");

  t=time(NULL);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("Calculation ended:%s\n", date);
  fprintf(parameters, "Calculation ended:%s\n", date);

  makePltFile(srcName, a, fileNamePrefix);
  printf("plt files are created.\n");

  free(bfst);
  free(blst);
  free(nxt);
  fprintf(stderr,"Bucket is freed\n");

  fclose(data);
  fclose(plot);
  fclose(partPlot);
  fclose(parameters);
  fclose(rigidBody);

  fclose(velocity);
  return  0;
  
}
