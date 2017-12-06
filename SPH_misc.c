#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"numbers.h"


void getSourceImageName(FILE *fp, char srcName[])//fp is supposed to be numbers
{
  char temp[64];
  int charCount=0;
  int i;

  fgets(temp, 64, fp);

  for(i=0; i<64; i++){
    if(temp[i]=='.' && temp[i+1]=='p' && temp[i+2]=='n' && temp[i+3]=='g'){
      break;
    }
    if(temp[i]=='/'){
      printf("skipped ");
      continue;
    }
    if(temp[i]=='.'){
      printf (". skipped ");
      continue;
    }
    srcName[charCount]=temp[i];
    charCount++;
  }

}

void makeDatFileName(char fName[], char type[], char srcName[], double angVel)
{
  char directory[128];

  sprintf(directory, "./Source_%s/", srcName);
  if(angVel>=0){
    sprintf(fName, "%sangvel%.2f_dt%.5f_%s_%s.dat", directory, angVel, dt, srcName, type);
  }else if(angVel<0){
    sprintf(fName, "%sangvelmin%.2f_dt%.5f_%s_%s.dat", directory, -angVel, dt, srcName, type);
  }

}

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
    fprintf(fp,"%.3f %.3fe\n", p[i].px, p[i].py);
  }
  fprintf(fp, "\n\n");

}

void printBoundaryParticles(Particle_State p[], FILE *fp)
{
  int i;
  for(i=FLP; i<FLP+BP; i++){ 
    fprintf(fp,"%d %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %0.12e %d \n",i, p[i].px, p[i].py,
            p[i].vx, p[i].vy, sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy),
            p[i].rho, p[i].p, p[i].ax, p[i].ay, p[i].color);
  }
  fprintf(fp,"\n\n");
}

void printBoundaryPositions(Particle_State p[], FILE *fp)
{
  int i;

  for(i=FLP; i<FLP+BP; i++){ 
    fprintf(fp,"%.3f %.3f\n", p[i].px, p[i].py);
  }
  fprintf(fp, "\n\n");
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
    fprintf(fp,"%.3fe %.3f\n", p[i].px, p[i].py);
  }
  fprintf(fp,"\n\n");
}

void printParticlesAroundObstacle(Particle_State p[], FILE *fp, double com[])
{
  int i;
  double xrange=3.0;
  double yrange=3.0;
  double gx=com[0];
  double gy=com[1];

  fprintf(fp, "%.3f %.3f\n", p[0].px, p[0].py);
  for(i=1; i<FLP; i++){
    if(p[i].px>=gx-xrange && p[i].px<=gx+xrange){
      if(p[i].py>=gy-yrange && p[i].py<=gy+yrange){
        fprintf(fp, "%.3f %.3f\n", p[i].px, p[i].py);
      }
    }
  }
  fprintf(fp, "\n\n");
  for(i=FLP+BP; i<N; i++){
    fprintf(fp, "%.3f %.3f\n", p[i].px, p[i].py);
  }
  fprintf(fp, "\n\n");
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

void getCalculationRegion(double range[], Particle_State p[])
{
  int i;
  double xmin, xmax, ymin, ymax;
  xmin=(XSIZE*interval)/2.0;
  xmax=(XSIZE*interval)/2.0;
  ymin=(YSIZE*interval)/2.0;

  for(i=FLP; i<FLP+BP; i++){
    if(p[i].px<xmin){
      xmin=p[i].px;
    }
    if(p[i].px>xmax){
      xmax=p[i].px;
    }
    if(p[i].py<ymin){
      ymin=p[i].py;
    }
  }

  range[0]=xmin;
  range[1]=xmax;
  range[2]=ymin;

}

void makePltFile(char *srcName, double angVel, Particle_State p[]){
  FILE *plt;
  FILE *plt2;
  FILE *partPlot;

  char directory[128];
  char fName[128];
  char line[256];
  double range[3];
  getCalculationRegion(range, p);
  sprintf(directory, "./Source_%s/rigid_anime.plt", srcName);

  plt=fopen(directory,"w");
  if(plt==NULL){
    printf("open plt error\n");
  }
  sprintf(directory, "./Source_%s/plot.plt", srcName);
  plt2=fopen(directory, "w");
  if(plt2==NULL){
    printf("open plt2 error\n");
  }
  sprintf(directory, "./Source_%s/partPlot.plt", srcName);
  partPlot=fopen(directory,"w");
  if(plt2==NULL){
    printf("open plt2 error\n");
  }
  
  if(angVel>=0){
    sprintf(fName, "angvel%.2f_dt%.5f_%s_plot.dat", angVel, dt, srcName);
  }else if(angVel<0){
    sprintf(fName, "angvelmin%.2f_dt%.5f_%s_plot.dat", -angVel, dt,srcName);
  }

  sprintf(line,"set xrange[%.1f:%.1f]\n",(range[0]-1), (range[1]+1));
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  sprintf(line,"set yrange[%.1f:%.1f]\n",(range[2]-1), (YSIZE+10)*interval);
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  sprintf(line,"set size ratio 1.0 1.0\n");
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  if((int)((DAMPTIME/100)*2-20)<0){ 
    sprintf(line,"do for[i=1:t:2]{\n");
  }else{
    sprintf(line, "do for[i=%d:t:2]{\n", (int)((DAMPTIME/100)*2-20+1));
  }
  fprintf(plt,"%s",line);
  sprintf(line,"print i\n");
  fprintf(plt,"%s",line);
  sprintf(line,"plot '%s' index 0 u 1:2 w p pt 35, '%s' index i u 1:2 w p pt 22, '%s' index i+1 u 1:2 w p pt 18\n", fName, fName, fName);
  fprintf(plt,"%s",line);
  sprintf(line,"}\n");
  fprintf(plt,"%s",line);

  sprintf(line, "plot '%s' index 0 u 1:2 w p pt 35, '%s' index (t*200) u 1:2 w p pt 18, '%s' index (t*200+1) u 1:2 w p pt 22\n", fName, fName, fName);
  fprintf(plt2, "%s", line);

  
  if(angVel>=0){
    sprintf(fName, "angvel%.2f_dt%.5f_%s_partPlot.dat", angVel, dt, srcName);
  }else if(angVel<0){
    sprintf(fName, "angvelmin%.2f_dt%.5f_%s_partPlot.dat", -angVel, dt,srcName);
  }

  sprintf(line,"set xrange[%.1f:%.1f]\n",(range[0]-1), (range[1]+1));
  fprintf(partPlot,"%s",line);
  sprintf(line,"set yrange[%.1f:%.1f]\n",(range[2]-1), (YSIZE+10)*interval);
  fprintf(partPlot,"%s",line);
  sprintf(line,"set size ratio %.3f 1.0\n",(double)((YSIZE+20)/(XSIZE+20)));
  fprintf(partPlot,"%s",line);
  if((int)((DAMPTIME/100)*2-20)<0){ 
    sprintf(line,"do for[i=1:t:2]{\n");
  }else{
    sprintf(line, "do for[i=%d:t:2]{\n", (int)((DAMPTIME/100)*2-20+1));
  }
  fprintf(partPlot,"%s",line);
  sprintf(line,"print i\n");
  fprintf(partPlot,"%s",line);
  sprintf(line,"plot '%s' index 0 u 1:2 w p pt 35, '%s' index i u 1:2 w p pt 22, '%s' index i+1 u 1:2 w p pt 18\n", fName, fName, fName);
  fprintf(partPlot,"%s",line);
  sprintf(line,"}\n");
  fprintf(partPlot,"%s",line);

  fclose(plt);
  fclose(plt2);
  fclose(partPlot);
}
