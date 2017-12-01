//setting for dambreaking
#include"SPH.h"
#include"numbers.h"
#include<stdio.h>
#include<math.h>


void initialization(Particle_State p[], int particleNumber)	//make all values of particles zero
{
  int	i;
  for(i=0; i<particleNumber; i++){
    p[i].inRegion = 0;//1 means particle is in region for calculatin 
    p[i].color=0;//1 means surface is hydrophilic, 2 means surface is hydrophilic
    p[i].px  = -100;
    p[i].py  = -100; 
    p[i].prepx=-100;
    p[i].prepy=-100;
    p[i].vx  = 0;
    p[i].vy  = 0;
    p[i].vxh = 0;
    p[i].vyh = 0;
    p[i].ax  = 0;
    p[i].ay  = 0;
    p[i].rho = 0;
    p[i].p   = 0;
    p[i].mass=m;
  }
}


int fluidParticles(Particle_State p[])//set fluid particles from fluid.txt
{
  FILE *fp;
  int i=0;
  int ret, px, py, r, gr, b;
  fp=fopen("fluid.txt","r");
  if(fp==NULL){
    printf("fluid.txt cannot be read.");
    return -1;
  }

  while((ret=fscanf(fp, "%d %d %*d %*d %*d", &px, &py))!=EOF){
      p[i].px=(px+1)*interval;
      p[i].py=(py+1)*0.088;//0.088 is used for exp2.png
      i++;
  }
  
  //fprintf(stderr,"%d\n", count);
  fclose(fp);
  return 0;
}


int wallParticles(Particle_State p[]){
  FILE *fp;
  int i=FLP;
  int ret, px,py;
  int green;
  fp=fopen("wall.txt","r");
  if(fp==NULL){
    printf("wall.txt cannot be read.");
    return -1;
  }
  while((ret = fscanf(fp, "%d %d %*d %d %*d", &px, &py, &green))!=EOF){
    p[i].px=(px+1)*interval;
    p[i].py=(py+1)*interval;
    if(green==255){
      p[i].color=1;//hydrophilic
    }else if(green==0){
      p[i].color=2;//hydrophobic
    }
    i++;
  }
  
  fclose(fp);
  return 0;
}


int obstacleBoundaryParticles(Particle_State p[])
{
  FILE *fp;
  int i=FLP+BP;
  int ret, px,py, green ;
  fp=fopen("obstacle.txt","r");
  if(fp==NULL){
    printf("obstacle.txt cannot be read.");
    return -1;
  }

  while((ret = fscanf(fp, "%d %d %*d %d %*d", &px, &py, &green))!=EOF){
    p[i].px=(px+1)*(interval);
    p[i].py=(py+1)*(interval);
    if(green==255){
      p[i].color=1;
    }else if(green==0){
      p[i].color=2;
    }
    i++;
  }
  
  fclose(fp);
  return 0;
}
