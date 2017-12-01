#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"numbers.h"



int main(int argc, char *argv[]){
    
  Particle_State a[N];
  int *bfst, *blst, *nxt;
  int countPer=0;
  int i;
  clock_t start, end;
  //  FILE *data;
  FILE *data;
  FILE *plot;
  FILE *parameters;
  FILE *numbers;
  FILE *rigidBody;
  char srcName[64];
  char fName[64];
  char date[64];
  time_t t=time(NULL);
  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("%s\n", date);
  double angVel=0;

  start=clock();  

  numbers=fopen("numbers.h","r");
  if(numbers==NULL){
    printf("Error!\n");
    return -1;
  }
  i=0;

  getSourceImageName(numbers, srcName);
  /*  fgets(temp, 100, numbers);
  int charcount=0;
  for(i=0; i<100; i++){
    if(temp[i]=='.' && temp[i+1]=='p' && temp[i+2]=='n' && temp[i+3]=='g'){
      break;
    }
    if(temp[i]=='/'){
      printf("skipped");
      continue;
    }
    if(temp[i]=='.'){
      printf(". skipped\n");
      continue;}
    srcName[charcount]=temp[i];
    charcount++;
    }*/

  printf("%s.png\n",srcName);
  if(argc==2){
    angVel=atof(argv[1]);
  }

  if(angVel>=0){
  sprintf(fName, "%f_%s_data.dat", angVel, srcName);
  }else if(angVel<0){
  sprintf(fName, "min%f_%s_data.dat", -angVel, srcName);
  }
  data=fopen(fName,"w");
  if(data==NULL){
    printf("File opening was failed.\n");
    return -1;
  }

  if(angVel>=0){
  sprintf(fName, "%f_%s_plot.dat", angVel, srcName);
  }else if(angVel<0){
  sprintf(fName, "min%f_%s_plot.dat", -angVel, srcName);
  }
  plot=fopen(fName,"w");
  if(plot==NULL){
    printf("File opening was failed.\n");
    return -1;
  }

  if(angVel>=0){
  sprintf(fName, "%f_%s_parameters.dat", angVel, srcName);
  }else if(angVel<0){
  sprintf(fName, "min%f_%s_parameters.dat", -angVel, srcName);
  }
  parameters=fopen(fName,"w");
  if(parameters==NULL){
    printf("File opening was failed.\n");
    return -1;
  }

  if(angVel>=0){
  sprintf(fName, "%f_%s_rigidBody.dat", angVel, srcName);
  }else if(angVel<0){
  sprintf(fName, "min%f_%s_rigidBody.dat", -angVel, srcName);
  }
  rigidBody=fopen(fName,"w");
  if(rigidBody==NULL){
    printf("File opening was failed. \n");
    return -1;
  }
  
  initialization(a, N);
  fprintf(stderr,"check initialization\n");
  fluidParticles(a);
  fprintf(stderr,"fluid Particles placed\n");
  wallParticles(a);
  fprintf(stderr,"Boundary Particles are placed\n");
  obstacleBoundaryParticles(a);
  fprintf(stderr, "Obstacke Particles are placed\n");
  allocateBucket(&bfst, &blst, &nxt);
  fprintf(stderr,"%d %d %d Bucket allocated\n", nBx, nBy, nBxy);
  checkParticle(a);
  
  
  makeBucket(bfst, blst, nxt, a);
  fprintf(stderr," makeBucket\n\n");

  fprintf(stderr,"Parameters:\nanglalrVelocity:%f\n", angVel);
  fprintf(stderr,"FLP=%d BP=%d OBP=%d\n", FLP, BP, OBP);
  fprintf(stderr,"m=%f h=%f rho0=%f dt=%f nu=%f g=%f gamm=%f T=%d DAMPTIME=%d\n\n\n",m,h,rho0,dt,nu,g,(double)gamm,T, DAMPTIME);

  fprintf(parameters,"Source Image %s.png\n",srcName);
  fprintf(parameters,"Anglar Velocity %f\n", angVel);
  fprintf(parameters,"FLP=%d BP=%d OBP=%d\n", FLP, BP, OBP);
  fprintf(parameters,"XSIZE=%d YSIZE=%d\n", XSIZE, YSIZE);
  fprintf(parameters,"m=%f h=%f rho0=%f dt=%f nu=%f g=%f gamm=%f T=%d DAMPTIME=%d\n\n\n",m,h,rho0,dt,nu,g,(double)gamm,T, DAMPTIME);
  
  rotateRigidBody(a, angVel);

  calcDensity(a, bfst, nxt);
  calcPressure(a);
  initializeAccel(a);
  calcAccelByExternalForces(a);
  calcAccelByPressure(a,bfst, nxt);
  calcAccelByViscosity(a,bfst, nxt,0);
  //  calcInterfacialForce(a, bfst, nxt);
  calcAccelByBoundaryForce(a, bfst, nxt);

  //calcAccelBySurfaceTension(a, bfst nxt);  
  //calcAccelByAdhesion(a, bfst, nxt);
  //printParticles(a,data);//Here shows parameters at t=0

  printBoundaryParticles(a, data);
  printFluidParticles(a, data);
  printObstacleParticles(a, data);

  printBoundaryPositions(a, plot);
  printFluidPositions(a, plot);
  printObstaclePositions(a,plot);


  //time development
  for(i=1; i<=T; i++){
    if(i==1){
      leapfrogStart(a);

    }else{
      leapfrogStep(a, i);
      }
    
    rigidBodyCorrection(a, rigidBody, i, angVel);
    checkParticle(a);
    makeBucket(bfst, blst, nxt, a);
    calcDensity(a, bfst, nxt);
    calcPressure(a);
    initializeAccel(a);
    calcAccelByExternalForces(a);
    calcAccelByPressure(a,bfst, nxt);
    calcAccelByViscosity(a,bfst, nxt,i);
    calcInterfacialForce(a, bfst, nxt);
    calcAccelByBoundaryForce(a, bfst, nxt);
    //calcAccelBySurfaceTension(a, bfst, nxt);
    //calcAccelByAdhesion(a, bfst, nxt);
    
    if(i%100==0){
      printFluidParticles(a, data);//here show paremeters at t=(i*dt)
      printObstacleParticles(a, data);
      
      //      if(i/100<=(T/100)/2){
      printFluidPositions(a,plot);
      printObstaclePositions(a,plot);
        //}else if(i/100 > (int) ((T/100)/2)){
        //	printFluidPositions(a,plot2);
        //	printObstaclePositions(a,plot2);
        //      }

    }
    percentage(i, &countPer);
  }

  fprintf(stderr,"\n");

  makePltFile(srcName, angVel);

  free(bfst);
  free(blst);
  free(nxt);

  fclose(data);
  fclose(plot);
  fclose(parameters);
  fclose(rigidBody);

  end=clock();
  fprintf(stderr,"Processor time: %fs\n", (double)(end-start)/CLOCKS_PER_SEC);

  strftime(date, sizeof(date), "%Y/%m/%d %a %H:%M:%S", localtime(&t));
  printf("%s\n", date);

  return  0;
  
}
