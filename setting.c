//setting for dambreaking
#include"SPH.h"
#include"Parameters.h"
#include"numbers.h"
#include<stdio.h>
#include<math.h>


void initialization(Particle_State p[], RigidPreValue rig[])	//make all values of particles zero
{
  int	i;
  for(i=0; i<N; i++){
    p[i].inRegion = 0;//1 means particle is in region for calculatin 
    p[i].color=0;//1 means surface is hydrophilic, 2 means surface is hydrophobic
    p[i].px  = -100;
    p[i].py  = -100; 
    p[i].vx  = 0;
    p[i].vy  = 0;
    p[i].vxh = 0;
    p[i].vyh = 0;
    p[i].ax  = 0;
    p[i].ay  = 0;
    p[i].rho = rho0;
    p[i].p   = 0;
    p[i].mass=m;
  }
  for(i=0; i<OBP; i++){
    rig[i].prepx=-100;
    rig[i].prepy=-100;
  }
}


int fluidParticles(Particle_State p[])//set fluid particles from fluid.txt
{
  FILE *fp;
  int i=0;
  int ret, px, py;
  fp=fopen("fluid.txt","r");
  if(fp==NULL){
    printf("fluid.txt cannot be read.");
    return -1;
  }

  while((ret=fscanf(fp, "%d %d %*d %*d %*d", &px, &py))!=EOF){
      p[i].px=(px+1)*interval;
      p[i].py=(py+1)*interval*0.9;//0.088 is used for exp2.png
      i++;
  }
  
  //fprintf(stderr,"%d\n", count);
  fclose(fp);
  return 0;
}

int fluidParticlesFromDat(Particle_State p[])
{
  FILE *fp;
  int i=0;
  int lineCount=0;
  int particleNumber=0;
  float px, py, vx, vy, ax, ay, pres, rho;
  char readline[256];
  px=0, py=0, vx=0, vy=0, ax=0, ay=0, pres=0, rho=0;
  fp=fopen("initialFluidState.dat", "r");
  if(fp==NULL){
    printf(" cannot be read!");
    return -1;
  }

  while(fgets(readline, 256, fp)!=NULL){
    if(lineCount>=1){
      sscanf(readline, "%d %f %f %f %f %f %f %f %f ",&particleNumber, &px, &py, &vx, &vy, &ax, &ay, &rho, &pres);
    p[particleNumber].px=px;
    p[particleNumber].py=py;
    p[particleNumber].vx=vx;
    p[particleNumber].vy=vy;
    p[particleNumber].ax=ay;
    p[particleNumber].ay=ay;
    p[particleNumber].p=pres;
    p[particleNumber].rho=rho;
    p[particleNumber].vxh=vx-ax*dt/2.0;
    p[particleNumber].vyh=vy-ay*dt/2.0;
    i++;
    }
    lineCount++;
  }
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
    } else if(green==200){
      p[i].color=3;
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
       p[i].color=1;//hydrophilic
    }else if(green==0){
       p[i].color=2;//hydrophobic
     }
     else if(green==200){
       p[i].color=3;
     }
    p[i].mass=rigidMass;
    i++;
  }
  
  fclose(fp);
  return 0;
}

void setInitialVelocity(Particle_State p[], double impactVel)
{
  int i;
  double angle=ANGLE_OF_INCIDENT;
  for(i=FLP+BP; i<N; i++){
    //    p[i].py-=7.0;
    p[i].vxh+=impactVel*cos(angle);
    p[i].vyh+=-impactVel*sin(angle);
  }

}
