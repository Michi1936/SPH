#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"numbers.h"


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


void printFluidPositions(Particle_State p[], FILE *fp)
{
  int i;

  for(i=0; i<FLP; i++){ 
    fprintf(fp,"%0.12e %0.12e\n", p[i].px, p[i].py);
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

void printObstaclePositions(Particle_State p[], FILE *fp)
{
  int i;
  for(i=FLP+BP; i<N; i++){ 
    fprintf(fp,"%0.12e %0.12e\n", p[i].px, p[i].py);
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

void makePltFile(char *srcName, double angVel){
  FILE *plt;
  char fName[40];
  char line[256];
  plt=fopen("rigid_anime.plt","w");
  sprintf(fName, "%f_%s_plot.dat",angVel, srcName);

  sprintf(line,"set xrange[-1:%f]\n",(XSIZE+10)*interval);
  fprintf(plt,"%s",line);
  sprintf(line,"set yrange[-1:%f]\n",(YSIZE+10)*interval);
  fprintf(plt,"%s",line);
  sprintf(line,"set size ratio %f 1.0\n",(double)((YSIZE+20)/(XSIZE+20)));
  fprintf(plt,"%s",line);
  sprintf(line,"do for[i=1:t:2]{\n");
  fprintf(plt,"%s",line);
  sprintf(line,"print i\n");
  fprintf(plt,"%s",line);
  sprintf(line,"plot '%s' index 0 u 2:3 w p lt 5, '%s' index i u 1:2 w p lt 7, '%s' index i+1 u 1:2 w p lt 10\n", fName, fName, fName);
  fprintf(plt,"%s",line);
  sprintf(line,"}\n");
  fprintf(plt,"%s",line);
  fclose(plt);
}
