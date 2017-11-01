#include"SPH.h"
#include"numbers.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>


void printParticles(Particle_State p[], FILE *fp);
void percentage(int time, int *countPer);
void printBoundaryParticles(Particle_State p[], FILE *fp);
void printFluidParticles(Particle_State p[], FILE *fp);
void printObstacleParticles(Particle_State p[], FILE *fp);
void tipPosition(Particle_State p[], int time, FILE *tip);

void printParticles(Particle_State p[], FILE *fp)
{
  int i;
  for(i=0; i<FLP; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay);
  }
  fprintf(fp, "\n\n");
  for(i=FLP; i<FLP+BP; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay);
  }
  for(i=FLP+BP; i<N; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay);
  }
  fprintf(fp,"\n\n");
}


void printFluidParticles(Particle_State p[], FILE *fp)
{
  int i;

  for(i=0; i<FLP; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay);
  }
  fprintf(fp, "\n\n");

}

void printBoundaryParticles(Particle_State p[], FILE *fp)
{
  int i;
  for(i=FLP; i<FLP+BP; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay);
  }
  fprintf(fp,"\n\n");
}


void printObstacleParticles(Particle_State p[], FILE *fp)
{
  int i;
  for(i=FLP+BP; i<N; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay);
  }
  fprintf(fp,"\n\n");
}

void percentage(int time, int *countPer)
{
  double per;
  per = time/(double)T;
  per=(int)(per*100);
  if(per==*countPer){
    fprintf(stderr,"%d ", (int)per);
    *countPer=*countPer+1;
  }
}

void tipPosition(Particle_State p[], int time, FILE *tip)
{
  double max=0;
  int i;
  for(i=0; i<FLP; i++){
    if(p[i].px>max){
      max=p[i].px;
    }
  }
  fprintf(tip, "%f %f %f %f\n", dt*time, max-0.3, (dt*time)*sqrt(2.0*g/4.0), ((max-0.3)/4.0));
}

int main(void){
    
  Particle_State a[N];
  int *bfst, *blst, *nxt;
  int countPer=0;
  int i;
  clock_t start, end;
  start=clock();
  FILE *data;
  FILE *plot;
  FILE *paramTxt;
  FILE *tip;
  
  data=fopen("sample.dat", "w");
  paramTxt=fopen("parameters.dat","w");
  tip=fopen("tip.dat", "w");
  plot=fopen("plot.dat","w");

  if(data==NULL){
    printf("File opening was failed.\n");
    return -1;
  }
  if(paramTxt==NULL){
    printf("File opening was failed.\n");
    return -1;
  }
  if(tip==NULL){
    printf("File opening was failed.\n");
    return -1;
  }
  if(plot==NULL){
    printf("File opening was failed.\n");
    return -1;
  }


  initialization(a, N);
  fprintf(stderr,"check initialization\n");
  fluidParticles(a);
  fprintf(stderr,"check initialConditions\n");
  wallParticles(a);
  fprintf(stderr,"check wallParitcles\n");
  obstacleBoundaryParticles(a);
  fprintf(stderr, "check obstacle boundary\n");
  allocateBucket(&bfst, &blst, &nxt);
  fprintf(stderr,"%d %d %d check allocateBucket\n", nBx, nBy, nBxy);
  checkParticle(a);
  
  
  makeBucket(bfst, blst, nxt, a);
  fprintf(stderr,"check makeBucket\n");
  
  fprintf(stderr,"FLP=%d BP=%d OBP=%d\n", FLP, BP, OBP);
  fprintf(stderr,"m=%f h=%f rho0=%f dt=%f nu=%f g=%f gamm=%f T=%d\n\n\n",m,h,rho0,dt,nu,g,(double)gamm,T);
    
  fprintf(paramTxt,"FLP=%d BP=%d OBP=%d\n", FLP, BP, OBP);
  fprintf(paramTxt,"m=%f h=%f rho0=%f dt=%f nu=%f g=%f gamm=%f T=%d\n\n\n",m,h,rho0,dt,nu,g,(double)gamm,T);
  
  //calculations for rho0, p0, a0  
  calcDensity(a, bfst, blst, nxt);
  calcPressure(a);
  initializeAccel(a);
  calcAccelByPressure(a,bfst,blst, nxt);
  calcAccelByViscosity(a,bfst,blst, nxt);
  calcAccelByExternalForces(a,bfst, blst, nxt);
  //calcAccelBySurfaceTension(a, bfst, blst, nxt);
  calcAccelByBoundaryForce(a, bfst, nxt);
  //printParticles(a,data);//Here shows parameters at t=0
  printBoundaryParticles(a, plot);
  printFluidParticles(a, plot);
  printObstacleParticles(a, plot);
  tipPosition(a, 0, tip);

  //time development
  for(i=1; i<=T; i++){
    if(i==1){
      leapfrogStart(a);
    }else{
      leapfrogStep(a);
      }
    //rigidBodyCorrection(a, bfst, nxt);
    checkParticle(a);
    makeBucket(bfst, blst, nxt, a);
    calcDensity(a, bfst, blst, nxt);
    calcPressure(a);
    initializeAccel(a);
    calcAccelByPressure(a,bfst,blst, nxt);
    calcAccelByViscosity(a,bfst,blst, nxt);
    calcAccelByExternalForces(a,bfst, blst, nxt);
    //calcAccelBySurfaceTension(a, bfst, blst, nxt);
    calcAccelByBoundaryForce(a, bfst, nxt);
    //printParticles(a,data);
    tipPosition(a, i, tip);
    if(i%100==0){
    printFluidParticles(a, plot);//here show paremeters at t=(i*dt)
    printObstacleParticles(a, plot);
    // fprintf(stderr,"%d printed\n", i);
    }
    percentage(i, &countPer);
  }

  fprintf(stderr,"\n");
  free(bfst);
  free(blst);
  free(nxt);
  fclose(data);
  fclose(paramTxt);
  fclose(tip);
  fclose(plot);
  end=clock();
  fprintf(stderr,"Processor time: %fs\n", (double)(end-start)/CLOCKS_PER_SEC);
  return  0;
  
}
