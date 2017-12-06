#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"numbers.h"



int main(int argc, char *argv[]){
    
  Particle_State a[N];
  double com[2];
  int *bfst, *blst, *nxt;
  int countPer=0;
  int i;
  FILE *data;
  FILE *plot;
  FILE *partPlot;
  FILE *parameters;
  FILE *numbers;
  FILE *rigidBody;
  FILE *newFIle;
  char srcName[128];
  char fName[128];
  char date[128];
  char type[128];
  char direct[128];
  double angVel=0;


  time_t t=time(NULL);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("\n\nCalculation started:%s\n", date);

  //open numbers.h
  numbers=fopen("numbers.h","r");
  if(numbers==NULL){
    printf("numbers.h cannot be opened!\n");
    return -1;
  }

  getSourceImageName(numbers, srcName);
  printf("%s.png\n",srcName);
  
  if(argc==2){
    angVel=atof(argv[1]);
  }

  sprintf(type, "plot");
  makeDatFileName(fName, type, srcName, angVel);
  if((newFIle=fopen(fName, "r"))!=NULL){
    int bool;
    fprintf(stderr, "\nWarning!\n");
    fprintf(stderr, "File %s is already exists.\nWould you like to overwrite this file?\n", fName);
    fprintf(stderr, "Yes:1   No:0\n");
    scanf("%d", &bool);
    if(bool==0){
      return -1;
    }
}

  plot=fopen(fName,"w");
  if(plot==NULL){
    printf("plot cannot be opened!\n");
    return -1;
  }

  sprintf(type, "data");
  makeDatFileName(fName, type, srcName, angVel);
  data=fopen(fName,"w");
  if(data==NULL){
    printf("data cannot be opened!\n");
    return -1;
  }

  sprintf(type, "parameters");
  makeDatFileName(fName, type, srcName, angVel);
  parameters=fopen(fName,"w");
  if(parameters==NULL){
    printf("parameters cannot be opened!\n");
    return -1;
  }

  sprintf(type, "rigidBody");
  makeDatFileName(fName, type, srcName, angVel);
  rigidBody=fopen(fName,"w");
  if(rigidBody==NULL){
    printf("rigidBody cannot be opened!\n");
    return -1;
  }

  sprintf(type, "partPlot");
  makeDatFileName(fName, type, srcName, angVel);
  partPlot=fopen(fName,"w");
  if(partPlot==NULL){
    printf("partPlot cannot be opened!\n");
    return -1;
  }
  
  //placing particles
  initialization(a, N);
  fprintf(stderr,"check initialization\n");
  fluidParticles(a);
  fprintf(stderr,"fluid Particles placed\n");
  wallParticles(a);
  fprintf(stderr,"Boundary Particles are placed\n");
  obstacleBoundaryParticles(a);
  fprintf(stderr, "Obstacke Particles are placed\n");


  //allocating bucket
  allocateBucket(&bfst, &blst, &nxt);
  fprintf(stderr,"%d %d %d Bucket allocated\n", nBx, nBy, nBxy);
  checkParticle(a);  
  makeBucket(bfst, blst, nxt, a);
  fprintf(stderr,"makeBucket\n\n");

  //printint parameters-------------------
  fprintf(stderr,"Parameters:\nanglalrVelocity:%f\n", angVel);
  fprintf(stderr,"FLP=%d BP=%d OBP=%d\n", FLP, BP, OBP);
  fprintf(stderr,"m=%f h=%f rho0=%f dt=%f nu=%f g=%f gamm=%f T=%d DAMPTIME=%d\n\n\n",m,h,rho0,dt,nu,g,(double)gamm,T, DAMPTIME);
  fprintf(parameters,"Source Image %s.png\n",srcName);
  fprintf(parameters,"Anglar Velocity %f\n", angVel);
  fprintf(parameters,"FLP=%d BP=%d OBP=%d\n", FLP, BP, OBP);
  fprintf(parameters,"XSIZE=%d YSIZE=%d\n", XSIZE, YSIZE);
  fprintf(parameters,"m=%f h=%f rho0=%f dt=%f nu=%f g=%f gamm=%f T=%d DAMPTIME=%d\n\n\n",m,h,rho0,dt,nu,g,(double)gamm,T, DAMPTIME);
  fprintf(parameters, "Calculation started%s\n", date);

  //calculating initial state
  rotateRigidBody(a, angVel);
  calcDensity(a, bfst, nxt);
  calcPressure(a);
  initializeAccel(a);
  calcAccelByExternalForces(a);
  calcAccelByPressure(a,bfst, nxt);
  calcAccelByViscosity(a,bfst, nxt,0);
  //  calcInterfacialForce(a, bfst, nxt);
  calcAccelByBoundaryForce(a, bfst, nxt);


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


  //time development
  for(i=1; i<=T; i++){
    if(i==1){
      leapfrogStart(a);
    }else{
      leapfrogStep(a, i);
      }

    rigidBodyCorrection(a, rigidBody, i, angVel, com);
    checkParticle(a);
    makeBucket(bfst, blst, nxt, a);

    calcDensity(a, bfst, nxt);
    calcPressure(a);
    initializeAccel(a);
    calcAccelByExternalForces(a);
    calcAccelByPressure(a,bfst, nxt);
    calcAccelByViscosity(a,bfst, nxt,i);
    //    calcInterfacialForce(a, bfst, nxt);
    calcAccelByBoundaryForce(a, bfst, nxt);
    
    if(i%100==0){
      printFluidParticles(a, data);//here show paremeters at t=(i*dt)
      printObstacleParticles(a, data);

      printFluidPositions(a,plot);
      printObstaclePositions(a,plot);
      
      printParticlesAroundObstacle(a, partPlot, com);
    }
    percentage(i, &countPer);
  }

  fprintf(stderr,"\n");

  makePltFile(srcName, angVel, a);

  free(bfst);
  free(blst);
  free(nxt);

  t=time(NULL);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("Calculation ended:%s\n", date);
  fprintf(parameters, "Calculation ended%s\n", date);

  fclose(data);
  fclose(plot);
  fclose(partPlot);
  fclose(parameters);
  fclose(rigidBody);

  return  0;
  
}
