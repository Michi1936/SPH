#include"SPH.h"
#include"numbers.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
//Equations below are based upon "Simulations of single bubbles rising through viscous liquids using SPH" by Szewc, Pozorski, Minier ,2013

//Monaghan(2005) cubic spline is used 
double cubicSpline1(double q)//cubic spline for 0<=q<=1
{
  return (15.0/(pow(h,2)*14.0*M_PI))*(pow(2.0-q,3)-4.0*pow(1.0-q,3));
  //  return (315.0/(64.0*M_PI*pow(h,9)))*pow(h*h-q*q,3);
}

double cubicSpline2(double q)//cubic spline for 1<=q<=2
{
  return (15.0/(pow(h,2)*14.0*M_PI))*(pow(2.0-q, 3));
}

double kernel(Particle_State p1, Particle_State p2)//p1 is central particle
{
  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  if(0<=q && q<=1){
    return cubicSpline1(q);
  }
   else if(1<=q && q<=2){
    return cubicSpline2(q);
    }
    else {
    return 0;
  }

}
double Wendrand(Particle_State p1, Particle_State p2){

  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double val=0;
  if(q<2){
    val = 21*pow(1.0-(q/2.0),4)*(2.0*q+1.0)/(16*M_PI*h*h*h);
  }else{
    val=0;
  }
  return val;
}

double gradKernel(Particle_State p1, Particle_State p2, int axis)//calculate gradient of kernel
{
  
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double coeff_x = (dx)/(dist*h+epsilon);
  double coeff_y = (dy)/(dist*h+epsilon);
  double dh = 0.01;
  if(axis==0){//x direction 
    if(0<=q && q<=1){
      return (15.0/(pow(h,2)*14.0*M_PI))*coeff_x*(12.0*pow(1.0-q,2)-3.0*pow(2.0-q,2));
    }
    else if(1<=q && q<=2){
      return -(15.0/(pow(h,2)*14.0*M_PI))*coeff_x*(3.0*pow(2.0-q,2));
    }
    else {
      return 0;
    }
  }
  
  else if(axis==1){//y direction 
    if(0<=q && q<=1){
      return (15.0/(pow(h,2)*14.0*M_PI))*coeff_y*(12.0*pow(1.0-q,2)-3.0*pow(2.0-q,2));
    }
    else if(1<=q && q<=2){
      return -(15.0/(pow(h,2)*14.0*M_PI))*coeff_y*(3.0*pow(2.0-q,2));
    }
    else {
      return 0;
    }
  }
  else{
    return 0;
  }
}

double gradWendlandKernel(Particle_State p1, Particle_State p2, int axis)
{
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double coeff_x = (dx)/(dist*h+epsilon);
  double coeff_y = (dy)/(dist*h+epsilon);
  double val=0;
  if(axis==0){//x direction 
    if(q<2){
      val = (21.0/(16.0*M_PI*pow(h,3)))*coeff_x*2.0*pow(1.0-(q/2.0),3)*(-5.0*q/2.0);
    }
    else{
      val=0;
    }
  }
    else if(axis==1){//y direction 
    if(q<2){
      val = (21.0/(16.0*M_PI*pow(h,3)))*coeff_y*2.0*pow(1.0-(q/2.0),3)*(-5.0*q/2.0);
    }
    else{
      val=0;
    }
  }
  else{
    return 0;
  }
  return val;
}

//Muller(2003)  Kernel for viscosity term is used
double Laplacian(Particle_State p1, Particle_State p2)
{
  double ans=0;
  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  if(0<= dist && dist <= h){
    ans=(45.0/(M_PI*pow(h,6)))*(h-dist);
  }
  else{
    ans = 0;
  }
  return ans;
}

double cohesion(double r)//caluculating cohesion term 
{
  double coef=32.0/(M_PI*pow(h,9));
  if(2*r>h && r<=h){
    return (coef*pow(h-r,3)*pow(r,3));
  }
  else if(r>0 && 2*r<=h){
    return (coef*(2.0*pow(h-r,3)*pow(r,3) - pow(h,6)/64.0));
  }
  else return 0;
}

void calcDensity(Particle_State p[], int bfst[], int blst[], int nxt[])
{
  int i;
  for(i=0; i<N; i++){
    p[i].rho=0;
  }

  //  fprintf(stderr,"density initialized\n");
  for(i=0; i<N; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      //fprintf(stderr, "%f %f %d %d %f",p[i].px ,p[i].py, ix, iy, BktLgth );
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb = 0;
          jb = jx + jy*nBx;
          int j = bfst[jb];
          //fprintf(stderr,"%d bfst accessed, %d %d %d\n", jb, i, jx, jy);
          while(j!=-1){
            double rhoij=0;
            rhoij=kernel(p[i], p[j]);
            p[i].rho+=p[i].mass*rhoij;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }
}


//Muller(2005) pressure model is used
void calcPressure(Particle_State p[])
{
  //set p=0
  int i;
  double coef;
  double cs=88.5;
  coef=(rho0*pow(cs,2))/7.0;
  for(i=0; i<N; i++){
    p[i].p=0;
  }
  for(i=0; i<N; i++){

    p[i].p = coef*(pow(p[i].rho/rho0,7)-1.0);//Tait equation
    if(p[i].p<0){
      p[i].p=0;
    }

  }

}
void initializeAccel(Particle_State p[])
{

  int i;
  for(i=0; i<N; i++){
    p[i].ax=0;
    p[i].ay=0;
  }
}

void calcInverseParticleVolume(Particle_State p[], int bfst[], int nxt[],double Theta[])
{
  int i,j;

  //inverse of the particle volume(Hu and Adams,2006)
  //  double Theta[FLP+BP];
  for (i=0; i<FLP+BP; i++){
    Theta[i]=0;
  }

  for(i=0; i<FLP+BP; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      //fprintf(stderr, "%f %f %d %d %f",p[i].px ,p[i].py, ix, iy, BktLgth );
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb = 0;
          jb = jx + jy*nBx;
          int j = bfst[jb];
          //fprintf(stderr, "j=%d, i=%d, ix=%d, iy=%d, jx=%d, jy=%d\n", j, i, ix, iy, jx, jy);
          while(j!=-1){
            Theta[i]+=kernel(p[i], p[j]);
            j = nxt[j];
            //      fprintf(stderr,"Theta_i is summed\n");
          }
        }
      }
    }else{Theta[i]=0;}
  }

}

void calcAccelByExternalForces(Particle_State p[], int bfst[], int blst[], int nxt[])
{
  int i;
for(i=0; i<FLP; i++){
    double aijx, aijy;
    aijx     = 0;
    aijy     = - g;//gravitational force
    p[i].ax += aijx;
    p[i].ay += aijy;
  }
}

void calcAccelByPressure(Particle_State p[], int bfst[], int blst[], int nxt[], double Theta[])
{
  int i;

  for(i=0; i<FLP+BP; i++){  
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb = jx + jy*nBx;
          int j = bfst[jb];
          if(j==-1)continue;
          for(;;){
            double aijx=0;
            double aijy=0;
            double presCoef=0;
           presCoef=(p[i].p/(pow(Theta[i],2)+epsilon)+p[j].p/(pow(Theta[j],2)+epsilon))/(p[i].mass+epsilon);

           aijx=presCoef*gradKernel(p[i], p[j], 0);
           aijy=presCoef*gradKernel(p[i], p[j], 1);
           p[i].ax += aijx;
           p[i].ay += aijy;
           j = nxt[j];
           if(j==-1)break;
          }
        }
      }
    } 
  }
}


void calcAccelByViscosity(Particle_State p[], int bfst[], int blst[], int nxt[], double Theta[])
    //Muller(2005) Weakly compressible SPH for free surface flow model is used.
{
  int i;
 
  for(i=0; i<FLP+BP; i++){

        if(p[i].inRegion==1){
          int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
          int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
          //fprintf(stderr, "%f %f %d %d %f",p[i].px ,p[i].py, ix, iy, BktLgth );
          int jx, jy;
          for(jx=ix-1; jx<=ix+1; jx++){
            for(jy=iy-1; jy<=iy+1; jy++){
              //fprintf(stderr,"%d bfst accessed, %d %d %d\n", jb, i, jx, jy);
              int jb=0;
              jb=jx+jy*nBx;
              int j = bfst[jb];
              double aijx, aijy;
              double viscCoef;
              aijx=0, aijy=0;
              viscCoef=0;
              double dx = fabs(p[i].px-p[j].px);
              double dy = fabs(p[i].py-p[j].py);
              double rx = (p[i].px-p[j].px);
              double ry = (p[i].py-p[j].py);
              double dist = dx*dx+dy*dy+epsilon;
          
          viscCoef=( (2.0*p[i].mu*p[j].mu)/(p[i].mu+p[j].mu+epsilon) )*(1.0/(pow(Theta[i],2)+epsilon)+1.0/(pow(Theta[j],2)+epsilon))*(rx*gradKernel(p[i],p[j],0)+ry*gradKernel(p[i],p[j],1))/(p[i].mass*(pow(dist,2)+pow(0.01,2)));
          aijx = viscCoef*(p[i].vx-p[i].vy);
          aijy = viscCoef*(p[i].vy-p[j].vy);
          p[i].ax+=aijx;
          p[i].ay+=aijy;
          j = nxt[j];
        }
      }
      //      fprintf(stderr, "%f  %f \n", p[i].ax, p[i].ay);
    }
  }
}

void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int blst[], int nxt[])
{
  int i,j;
  
  double *nx, *ny;
    nx=(double*)malloc(sizeof(double)*FLP);
    ny=(double*)malloc(sizeof(double)*FLP);
    int k;
    for(i=0; i<FLP; i++){
      nx[i]=0;
      ny[i]=0;
    }
    for(i=0; i<FLP; i++){
      for(k=0; k<FLP; k++){
        nx[i]+=h*m*gradKernel(p[i],p[k],0)/(p[k].rho+epsilon);
        ny[i]+=h*m*gradKernel(p[i],p[k],1)/(p[k].rho+epsilon);
      }
    }

    //calculation of cohesion force and surface area minimization
    //Akinci(2013) cohesion and surface tension model is used 
    for(i=0; i<FLP; i++){
      if(p[i].inRegion==1){
        int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
        int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
        //fprintf(stderr, "%f %f %d %d %f",p[i].px ,p[i].py, ix, iy, BktLgth );
        int jx, jy;
        
        for(jx=ix-1; jx<=ix+1; jx++){
          for(jy=iy-1; jy<=iy+1; jy++){
            int jb = jx + jy*nBx;
            int j = bfst[jb];
            //fprintf(stderr,"%d bfst accessed, %d %d %d\n", jb, i, jx, jy);
            if(j==-1)continue;
            for(;;){
              double aijx, aijy;
              double dx, dy, dist, Kij;
              
              aijx=0, aijy=0;
              dx = fabs(p[i].px-p[j].px);
              dy = fabs(p[i].py-p[j].py);
              
              dist = sqrt(dx*dx+dy*dy)+epsilon;
              Kij=2.0*rho0/(p[i].rho+p[j].rho+epsilon);
              
              aijx = -gamm*m*(m*cohesion(dist)*dx/dist + (nx[i]-nx[j]));
              aijy = -gamm*m*(m*cohesion(dist)*dy/dist + (ny[i]-ny[j]));
              
              aijx=aijx/(p[i].rho+epsilon);
              aijy=aijy/(p[i].rho+epsilon);
              p[i].ax+=aijx;
              p[i].ay+=aijy;
              
              j = nxt[j];
              if(j==-1)break;
            }
          }
        }
      }
    }
    // fprintf(stderr,"Acceleration is calculated\n");
    free(nx);
    free(ny);
  
}

void timeDevelopment(Particle_State p[])
{
  int i;
  for(i=0; i<FLP; i++){
    //Symplectic Euler method
    p[i].vx += p[i].ax*dt;
    p[i].vy += p[i].ay*dt;
    p[i].px+=p[i].vx*dt; 
    p[i].py+=p[i].vy*dt; 
  }
  //  fprintf(stderr, "time developed\n")
}


void leapfrogStart(Particle_State p[])
{
  int i;
  for(i=0; i<FLP; i++){
 
    p[i].vxh=p[i].vx+p[i].ax*dt/2.0;
    p[i].vyh=p[i].vy+p[i].ay*dt/2.0;
    p[i].vx+=p[i].ax*dt;
    p[i].vy+=p[i].ay*dt;
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
    if(i==2){

      //     fprintf(stderr, "%0.12e %0.12e %0.12e %0.12e %0.12e %0.12e\n", p[i].ax, p[i].ay, p[i].vxh, p[i].vyh, p[i].px, p[i].py);
    }
  }
}

void boundaryCondition(Particle_State p[])
{
  double LW=interval*3;
  double RW=4.0-interval*2;
  double BW=0.05*3;
  int i;
  for(i=0; i<FLP; i++){
    if(p[i].px<LW || p[i].px>RW){
      p[i].vx=0;
      p[i].vy=0;
    }

    if(p[i].py<BW){
      p[i].vx=0;
      p[i].vy=0;
    }
  }
          
}

void leapfrogStep(Particle_State p[])
{
  int i;
  for(i=0; i<FLP; i++){
    p[i].vxh+=p[i].ax*dt;
    p[i].vyh+=p[i].ay*dt;
    p[i].vx=p[i].vxh+p[i].ax*dt/2.0;
    p[i].vy=p[i].vyh+p[i].ay*dt/2.0;
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
  }
  
}


