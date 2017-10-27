#include"SPH.h"
#include"numbers.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>


void printParticles(Particle_State p[], FILE *fp);
void percentage(int time, int *countPer);


    
void printParticles(Particle_State p[], FILE *fp){
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

  //  printf("Printed\n");
}

void percentage(int time, int *countPer)
{
  double per;
  per=(int)(per*100);
  if(per==*countPer){
    fprintf(stderr,"%d ", (int)per);
    *countPer=*countPer+1;
  }
 
  
}

int main(void){
    
  Particle_State a[N];
  double Theta[FLP+BP];
  int *bfst, *blst, *nxt;
  int countPer=0;
  int i;
  clock_t start, end;
  start=clock();
  FILE *fp;
  FILE *paramTxt;
  fp=fopen("sample.dat", "w");
  paramTxt=fopen("parameters.dat","w");

  if(fp==NULL){
    printf("File opening was failed.\n");
    return -1;
  }

  if(paramTxt==NULL){
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
  fprintf(stderr,"check checkParticle");

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
  calcInverseParticleVolume(a, bfst, nxt, Theta);
  calcAccelByPressure(a,bfst,blst, nxt,Theta);
  calcAccelByViscosity(a,bfst,blst, nxt,Theta);
  calcAccelByExternalForces(a,bfst, blst, nxt);
  calcAccelBySurfaceTension(a, bfst, blst, nxt);
  
  printParticles(a,fp);//Here shows parameters at t=0
  
  //time development
  for(i=1; i<=T; i++){
    if(i==1){
      leapfrogStart(a);
    }else{
      leapfrogStep(a);
      //      boundaryCondition(a);
    }
    checkParticle(a);
    makeBucket(bfst, blst, nxt, a);
    calcDensity(a, bfst, blst, nxt);
    calcPressure(a);
    initializeAccel(a);
    calcInverseParticleVolume(a, bfst, nxt, Theta);
    calcAccelByPressure(a,bfst,blst, nxt, Theta);
    calcAccelByViscosity(a,bfst,blst, nxt, Theta);
    calcAccelByExternalForces(a,bfst, blst, nxt);
    calcAccelBySurfaceTension(a, bfst, blst, nxt);
    printParticles(a, fp);//here show paremeters at t=(i*dt)
    percentage(i, &countPer);
  }

  fprintf(stderr,"\n");
  free(bfst);
  free(blst);
  free(nxt);
  fclose(fp);
  end=clock();
  fprintf(stderr,"Processor time: %fs\n", (double)(end-start)/CLOCKS_PER_SEC);
  return  0;
  
}
