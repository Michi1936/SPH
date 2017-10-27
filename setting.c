//setting for dambreaking
#include"SPH.h"
#include"numbers.h"
#include<stdio.h>
#include<math.h>


void initialization(Particle_State p[], int particleNumber)	//make all values of particles zero
{
  int	i;
  int half=particleNumber/2;

  for(i=0; i<particleNumber; i++){
    p[i].inRegion = 0;//1 means particle is in region for calculatin 
    p[i].px  = -100;
    p[i].py  = -100; 
    p[i].vx  = 0;
    p[i].vy  = 0;
    p[i].vxh = 0;
    p[i].vyh = 0;
    p[i].ax  = 0;
    p[i].ay  = 0;
    p[i].mass= 1.0;
    p[i].rho = 0;
    p[i].p   = 0;
    p[i].mu=0.3;
  }
}


int fluidParticles(Particle_State p[])//set fluid particles from fluid.txt
{
  FILE *fp;
  int i=0;
  int ret, count=0, px, py, r, gr, b;
  fp=fopen("fluid.txt","r");
  if(fp==NULL){
    printf("fluid.txt cannot be read.");
    return -1;
  }

  while((ret=fscanf(fp, "%d %d %d %d %d", &px, &py, &r, &gr, &b))!=EOF){
      p[i].px=(px+1)*interval;
      p[i].py=(py+1)*interval;
      p[i].vy=(gr-127)*0.2;
      i++;
  }
  /*
    fp=fopen("fluid.txt","r");
    for(i=0; i<FLP; i++){
      fscanf(fp,"%d %d %d %d %d", &px, &py, &r, &gr, &b);
      p[i].px=(px+1)*interval;
      p[i].py=(py+1)*interval;
      p[i].vy=(gr-127)*0.2;
      fprintf(stderr, "vy=%f %d", p[i].vy, gr);
      //      fprintf(stderr,"%f %f %f\n", p[i].px, p[i].py, p[i].vy);
      */
  
  //fprintf(stderr,"%d\n", count);
  fclose(fp);
  return 0;
}


int wallParticles(Particle_State p[]){
  FILE *fp;
  int i=FLP;
  int ret, count=0, px,py;
  fp=fopen("wall.txt","r");
  if(fp==NULL){
    printf("wall.txt cannot be read.");
    return -1;
  }
  while((ret = fscanf(fp, "%d %d %*d %*d %*d", &px, &py))!=EOF){
    p[i].px=(px+1)*interval;
    p[i].py=(py+1)*interval;
    i++;
  }
  
  /*  fp=fopen("wall.txt","r");
  for(i=FLP; i<FLP+BP; i++){
    fscanf(fp,"%d %d %d %d %d", &px, &py);
    p[i].px=(px+1)*interval;
    p[i].py=(py+1)*interval;
    fprintf(stderr,"%f %f\n", p[i].px, p[i].py);
    
    }*/
  fprintf(stderr,"%d\n", count);
  fclose(fp);
  return 0;
}


void obstacleBoundaryParticles(Particle_State obp[])
{
  int i;
  double dif = interval;
  double x,y;
  x=1.0-dif;
  y=0.25;
  int count=0;
  /*  for(i=0; i<OBP; i++){
      obp[i].px=x;
      obp[i].py=y;
      count++;
      x+=dif;
      if(count==3){
        y+=dif;
        x=1.0-dif;
        count=0;
      }
      }*/
 
      for(i=FLP+BP; i<N; i++){
      obp[i].px=-100;
      obp[i].py=-100;
      }
}
