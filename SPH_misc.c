#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <string.h>
#include"SPH.h"
#include"numbers.h"


void getSourceImageName(FILE *fp, char srcName[])//fp is supposed to be numbers
    //extracting source image name on 1st line of numbers.h
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
    srcName[charCount]=temp[i];
    fprintf(stderr ,"%c ", srcName[charCount]);
    charCount++;
  }
  fprintf(stderr, "%d\n", charCount);
  srcName[charCount]='\0';
  fprintf(stderr, "\nsrcName is %s\n", srcName);
  printf("\n");
}

double calcRadius(Particle_State p[])
{
  double gx, gy;
  double qx[OBP], qy[OBP];
  double radius;
  int i;
  gx=0, gy=0;

  for(i=0; i<OBP; i++){
    qx[i]=0;
    qy[i]=0;
  }
  
  for (i=FLP+BP; i<N; i++){//calculating center of mass
    gx+=p[i].px/OBP;
    gy+=p[i].py/OBP;
  }
  
  for (i=FLP+BP; i<N; i++){//calculating vector between center of mass and ith particle
    qx[i-FLP-BP]=p[i].px-gx;
    qy[i-FLP-BP]=p[i].py-gy;
  }

  for(i=FLP+BP; i<N; i++){//calculating inertia
    double dist=sqrt(qx[i-FLP-BP]*qx[i-FLP-BP] + qy[i-FLP-BP]*qy[i-FLP-BP]);
    if(dist>radius){
      radius=dist;
    }
  }
  fprintf(stderr, "radius=%f\n", radius);
  return radius;

}

void makeFileNamePrefix(char fileNamePrefix[], char srcName[], double angVel, double spinParam)
{
  char prefix[512];
  char tempFileName[512];
  char suffix[128];
  FILE *dammyPLOT;

  if(angVel>=0){
    sprintf(prefix, "angvel%.2f_dt%.6f_sp%.3f_nu%.5f_%s",  angVel, dt, spinParam,nu, srcName);
  }else if(angVel<0){
    sprintf(prefix, "angvelmin%.2f_dt%.6f_sp%.3f_nu%.5f_%s",  -angVel, dt, spinParam,nu,srcName);
  }

  sprintf(tempFileName,"./Source_%s/%s_plot.dat", srcName,prefix);
  sprintf(fileNamePrefix, "%s", prefix);
  if((dammyPLOT=fopen(tempFileName,"r"))!=NULL){
    int bool;
    fprintf(stderr, "\nWarning!\n");
    fprintf(stderr, "File %s already exists.\n", tempFileName);
    fprintf(stderr,"1:Overwrite this file\n");
    fprintf(stderr,"2:Rename this file\n");
    fprintf(stderr, "0:Stop this programm\n");
    scanf("%d", &bool);
    if(bool==1){
      sprintf(fileNamePrefix, "%s", prefix);
      fprintf(stderr, "bool=%d\n", bool);
    }else if(bool==2){
      fprintf(stderr,"Enter suffix\n");
      scanf("%s", suffix);
      sprintf(fileNamePrefix, "%s_%s", prefix, suffix);
    }else{
      fprintf(stderr, "Program ended\n");
      exit(EXIT_FAILURE);
    }
  }

  fprintf(stderr, "makeFileName ended\n");
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

void printParameters(FILE *fp, double angVel, char srcName[], char date[], double spinParam)
{
  fprintf(stderr, "Parameters:\n");
  fprintf(stderr,"Source Image %s.png\n",srcName);
  fprintf(stderr,"Anglar Velocity:%f Impact Velocity:%f Angle of Incident:%f Spin Parameter:%f\n", angVel, IMPACT_VELOCITY, ANGLE_OF_INCIDENT, spinParam);
  fprintf(stderr,"FLP:%d BP:%d OBP:%d\n", FLP, BP, OBP);
  fprintf(stderr,"XSIZE:%d YSIZE:%d\n", XSIZE, YSIZE);
  fprintf(stderr, "FLUID_INTERACTION:%f HPHILY_INTERACTION:%F HPHOBY_INTERACTION:%F\n", FLUID_INTERACTION, HPHILY_INTERACTION, HPHOBY_INTERACTION);
  fprintf(stderr, "BOUNDARY_FORCE:%d\n", BOUNDARY_FORCE);
  fprintf(stderr,"m:%f rigidMass:%f h:%f rho0:%f dt:%f nu:%f g:%f T:%d \nDAMPTIME:%d MOTION_START_TIME:%d\n\n\n",m, rigidMass, h,rho0,dt,nu,g,T, DAMPTIME, MOTION_START_TIME);
  fprintf(stderr, "Calculation started:%s\n", date);
  
  fprintf(fp,"Source Image %s.png\n",srcName);
  fprintf(fp,"Anglar Velocity:%f Impact Velocity:%f Angle of Incident:%f Spin Parameter:%f\n", angVel, IMPACT_VELOCITY, ANGLE_OF_INCIDENT, spinParam);
  fprintf(fp,"FLP:%d BP:%d OBP:%d\n", FLP, BP, OBP);
  fprintf(fp,"XSIZE:%d YSIZE:%d\n", XSIZE, YSIZE);
  fprintf(fp, "FLUID_INTERACTION:%f HPHILY_INTERACTION:%F HPHOBY_INTERACTION:%F\n", FLUID_INTERACTION, HPHILY_INTERACTION, HPHOBY_INTERACTION);
  fprintf(fp, "BOUNDARY_FORCE:%d\n", BOUNDARY_FORCE);
  fprintf(fp,"m:%f rigidMass:%f h:%f rho0:%f dt:%f nu:%f g:%f T:%d DAMPTIME:%d MOTION_START_TIME:%d\n\n\n",m,rigidMass,h,rho0,dt,nu,g,T, DAMPTIME, MOTION_START_TIME);

  fprintf(fp, "Calculation started:%s\n", date);
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


void getMaxVelocity(Particle_State p[], FILE *fp, int time)
{
  int i;
  double max=0;
  for(i=0; i<FLP; i++){
    double vel=sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy);
    if(vel>max){
      max=vel;
    }
  }
  fprintf(fp, "%d %f\n", time, max);
}

void getCalculationRegion(double range[], Particle_State p[])
{
  int i;
  double xmin, xmax, ymin;
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

void makePltFile(char *srcName, Particle_State p[], char *fileNamePrefix)
{
  FILE *plt;
  FILE *plt2;
  FILE *partPlot;

  char directory[128];
  char fName[128];
  char line[1024];
  double range[3];
  double ratio=0;
  int sec;
  getCalculationRegion(range, p);
  ratio=(((YSIZE+10)*interval)-(range[2]-1))/(range[1]-range[0]+2);
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
  if(partPlot==NULL){
    printf("open plt2 error\n");
  }
  
  sprintf(fName, "%s_plot.dat", fileNamePrefix);
  
  sprintf(line,"set xrange[%.1f:%.1f]\n",(range[0]-1), (range[1]+1));
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  sprintf(line,"set yrange[%.1f:%.1f]\n",(range[2]-1), (YSIZE+10)*interval);
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  sprintf(line,"set size ratio %f\n",ratio);
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);

  //writing in rigid_anime.plt
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
  sprintf(line,"print %.2f\n", ratio);
  fprintf(plt,"%s",line);

  //writing in plot.plt
  sprintf(line, "plot '%s' index 0 u 1:2 w p pt 35, '%s' index (t*200) u 1:2 w p pt 22, '%s' index (t*200+1) u 1:2 w p pt 18\n", fName, fName, fName);
  fprintf(plt2, "%s", line);
  sprintf(line,"print %.2f\n", ratio);
  fprintf(plt2,"%s",line);


  //writing in partPlot.plt
  sprintf(fName, "%s_partPlot.dat", fileNamePrefix);
  sprintf(line,"set xrange[%.1f:%.1f]\n",(range[0]-1), (range[1]+1));
  fprintf(partPlot,"%s",line);
  sprintf(line,"set yrange[%.1f:%.1f]\n",(range[2]-1), (YSIZE+10)*interval);
  fprintf(partPlot,"%s",line);
  sprintf(line,"set size ratio %f\n",ratio);
  fprintf(partPlot,"%s",line);
  if((int)((DAMPTIME/100)*2-20)<0){ 
    sprintf(line,"do for[i=1:t:2]{\n");
  }else{
    sprintf(line, "do for[i=%d:t:2]{\n", (int)((DAMPTIME/100)*2-20+1));
  }
  fprintf(partPlot,"%s",line);
  sprintf(line,"print i\n");
  fprintf(partPlot,"%s",line);
  sprintf(line,"plot '%s' index 0 u 1:2 w p pt 35, '%s' index i u 1:2 w p pt 22, '%s' index i+1 u 1:2 w p pt 18\n", fName, fName,fName);
  fprintf(partPlot,"%s",line);
  sprintf(line,"}\n");
  fprintf(partPlot,"%s",line);
  sprintf(line,"print %.2f\n", ratio);
  fprintf(partPlot,"%s",line);

  fclose(plt);
  fclose(plt2);
  fclose(partPlot);
  fprintf(stderr, "plot files are made\n");
  

  //from here plt files are make in current derectory
  sprintf(directory, "./rigid_anime.plt");
  sprintf(fName, "./Source_%s/%s_plot.dat", srcName, fileNamePrefix);
  plt=fopen(directory,"w");
  if(plt==NULL){
    printf("open plt error\n");
  }

  sprintf(directory, "./plot.plt");
  plt2=fopen(directory, "w");
  if(plt2==NULL){
    printf("open plt2 error\n");
  }

  sprintf(directory, "./partPlot.plt");
  partPlot=fopen(directory,"w");
  if(partPlot==NULL){
    printf("open plt2 error\n");
  }
  
  sprintf(line,"set xrange[%.1f:%.1f]\n",(range[0]-1), (range[1]+1));
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  sprintf(line,"set yrange[%.1f:%.1f]\n",(range[2]-1), (YSIZE+10)*interval);
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);
  sprintf(line,"set size ratio %f\n",ratio);
  fprintf(plt,"%s",line);
  fprintf(plt2,"%s",line);

  //writing in rigid_anime.plt
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
  sprintf(line,"print %.2f\n", ratio);
  fprintf(plt,"%s",line);

  //writing in plot.plt
  sprintf(line, "plot '%s' index 0 u 1:2 w p pt 35, '%s' index (t*200) u 1:2 w p pt 22, '%s' index (t*200+1) u 1:2 w p pt 18\n", fName, fName, fName);
  fprintf(plt2, "%s", line);
  sprintf(line,"print %.2f\n", ratio);
  fprintf(plt2,"%s",line);


  //writing in partPlot.plt
  sprintf(fName, "./Source_%s/%s_partPlot.dat", srcName, fileNamePrefix);
  sprintf(line,"set xrange[%.1f:%.1f]\n",(range[0]-1), (range[1]+1));
  fprintf(partPlot,"%s",line);
  sprintf(line,"set yrange[%.1f:%.1f]\n",(range[2]-1), (YSIZE+10)*interval);
  fprintf(partPlot,"%s",line);
  if((int)((DAMPTIME/100)*2-20)<0){ 
    sprintf(line,"do for[i=1:t:2]{\n");
  }else{
    sprintf(line, "do for[i=%d:t:2]{\n", (int)((DAMPTIME/100)*2-20+1));
  }
  fprintf(partPlot,"%s",line);
  sprintf(line,"print i\n");
  fprintf(partPlot,"%s",line);
  sprintf(line,"plot '%s' index 0 u 1:2 w p pt 35, '%s' index i u 1:2 w p pt 22, '%s' index i+1 u 1:2 w p pt 18\n", fName, fName,fName);
  fprintf(partPlot,"%s",line);
  sprintf(line,"}\n");
  fprintf(partPlot,"%s",line);
  sprintf(line,"print %.2f\n", ratio);
  fprintf(partPlot,"%s",line);

  fclose(plt);
  fclose(plt2);
  fclose(partPlot);
  fprintf(stderr, "plot files are made in current derectory\n");
}
