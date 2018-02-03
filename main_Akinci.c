#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"Parameters.h"
#include"numbers.h"



int main(int argc, char *argv[]){
    
  Particle_State a[N];
  double com[2];
  double Psi[BP+OBP];
  int boundaryType[BP+OBP];
  RigidBodyValues rigV;
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
  char type[32];
  char fileNamePrefix[256];
  char date[256];
  double angVel=0;

  time_t t=time(NULL);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("\n\nCalculation started:%s\n", date);

  //defining anglar velocity
  if(argc==2){
    angVel=atof(argv[1]);
  }

  rigV.omega=angVel;

  //placing particles
  initialization(a, rigV);
  fprintf(stderr,"check initialization\n");
  fluidParticles(a);
  fprintf(stderr,"fluid Particles placed\n");
  wallParticles(a);
  fprintf(stderr,"Boundary Particles are placed\n");
  obstacleBoundaryParticles(a);
  fprintf(stderr, "Obstacke Particles are placed\n");

  //allocating bucket
  allocateBucket(&bfst, &blst, &nxt);
  fprintf(stderr,"%f %f %d %d %d Bucket allocated\n", BktLgth, BktNum, nBx, nBy, nBxy);
  checkParticle(a);  
  makeBucket(bfst, blst, nxt, a);
  fprintf(stderr,"makeBucket\n\n");

  spinParam=fabs(calcRadius(a)*angVel/(IMPACT_VELOCITY+epsilon));
  fprintf(stderr, "spin parameter:%f\n", spinParam);
  //open numbers.h
  numbers=fopen("numbers.h","r");
  if(numbers==NULL){
    printf("numbers.h cannot be opened!\n");
    return -1;
  }

  getSourceImageName(numbers, srcName);
  printf("%s.png\n",srcName);
  
  makeFileNamePrefix(fileNamePrefix, srcName, angVel, spinParam);
  fprintf(stderr, "\n\nprefix and suffix test\n./Source_%s/%s_plot.dat\n", srcName, fileNamePrefix);

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
  calcPsi(a, Psi, bfst, nxt, boundaryType);
  AkinciCalcDensity(a, Psi, bfst, nxt);
  calcPressure(a);

  initializeAccel(a);
  calcAccelByExternalForces(a);
  AkinciCalcAccelByPressure(a, Psi, bfst, nxt);
  AkinciCalcAccelByViscosity(a, Psi, bfst, nxt,0);

  if(FLUID_INTERACTION>epsilon){
    calcInterfacialForce(a, bfst, nxt);
  }
  if(BOUNDARY_FORCE==1){
    calcAccelByBoundaryForce(a, bfst, nxt);
  }
  

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
    rotateRigidBody(a, rigV, angVel);
  }
  
  //time development
  for(i=1; i<=T; i++){
    if(i==1){
      EulerCromerTimeIntegration(a);
    }else{
      EulerCromerTimeIntegration(a);
      if(i==MOTION_START_TIME){
        fprintf(stderr, "At %d rigid body is rotated\n", i);
        rotateRigidBody(a, rigV, angVel);
        setInitialVelocity(a);
        fprintf(stderr, "Rigid Body Velocity is set\n");
      }
    }

    rigidBodyTimeIntegration(a, rigV, rigidBody, i);
    checkParticle(a);
    makeBucket(bfst, blst, nxt, a);

    calcPsi(a, Psi, bfst, nxt, boundaryType);
    AkinciCalcDensity(a, Psi, bfst, nxt);
    calcPressure(a);

    initializeAccel(a);
    calcAccelByExternalForces(a);
    AkinciCalcAccelByPressure(a, Psi, bfst, nxt);
    AkinciCalcAccelByViscosity(a, Psi, bfst, nxt, i);
        
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
