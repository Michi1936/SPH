#include"SPH.h"
#include"numbers.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

double cubicSpline1(double q)//cubic spline for 0<=q<=1
{
  return Ch*(pow(2.0-q,3) - 4.0*pow(1.0-q,3));
}

double cubicSpline2(double q)//cubic spline for 1<=q<=2
{
  return Ch*pow(2.0-q,3);
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

double gradKernel(Particle_State p1, Particle_State p2, int axis)//calculate gradient of kernel
{
  
  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double coeff_x = (p1.px-p2.px)/(dist*h+epsilon);
  double coeff_y = (p1.py-p2.py)/(dist*h+epsilon);
  
  if(axis==0){//x direction 
    if(0<=q && q<=1){
      return coeff_x*((cubicSpline1(q+dh)-cubicSpline1(q-dh))/(2*dh));
    }
    else if(1<=q && q<=2){
      return coeff_x*((cubicSpline2(q+dh)-cubicSpline2(q-dh))/(2*dh));
    }
    else {
      return 0;
    }
  }
  else if(axis==1){//y direction 
    if(0<=q && q<=1){
      return coeff_y*((cubicSpline1(q+dh)-cubicSpline1(q-dh))/(2*dh));
    }
    else if(1<=q && q<=2){
      return coeff_y*((cubicSpline2(q+dh)-cubicSpline2(q-dh))/(2*dh));
    }
    else {
      return 0;
    }
  }
  else{
    return 0;
  }
}
  
double Laplacian(Particle_State p1, Particle_State p2)
{
  double ans=0;
  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  if(0<= dist && dist <= h){
    ans = 20/(3*M_PI*pow(h,5));
    ans = ans*(h-dist);
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
          if(j==-1)continue;
          for(;;){
            double rhoij=0;
            rhoij=m*kernel(p[i], p[j]);
            p[i].rho+=rhoij;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }
  //  fprintf(stderr, "Density calculated\n");
}


void calcPressure(Particle_State p[])
{
  //set p=0
  int i;
  for(i=0; i<N; i++){
    p[i].p=0;
  }
  
  for(i=0; i<N; i++){
    p[i].p = k1*(pow(p[i].rho/rho0,7)-1.0);//Tait equation
    if(p[i].p<0) {
      p[i].p=0;
    }
  }
}

void calcAcceleration(Particle_State p[], int bfst[], int blst[], int nxt[])
{
  int i,j;
  
  for(i=0; i<N; i++){
    p[i].ax=0;
    p[i].ay=0;
  }
 
  //calculation about external force
  for(i=0; i<FLP; i++){
    double aijx, aijy;
    aijx     = 0;
    aijy     = - g;//gravitational force
    p[i].ax += aijx;
    p[i].ay += aijy;
  }

  //calculation about viscosity
  //Muller(2003) viscosity model is used 
  for(i=0; i<FLP+BP; i++){
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
            aijx=0, aijy=0;
            double viscCoef=0;
            double dx = fabs(p[i].px-p[j].px);
            double dy = fabs(p[i].py-p[j].py);
            double dist = dx*dx+dy*dy+epsilon;
            viscCoef = nu*m*Laplacian(p[i], p[j])/p[j].rho;
            aijx = -viscCoef*(p[i].vx-p[j].vx);
            aijy = -viscCoef*(p[i].vy-p[j].vy);
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }
  //calculation about pressure   
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
            double aijx=0, aijy=0;
            aijx = -m*(p[i].p/(pow(p[i].rho,2)+epsilon) + p[j].p/(pow(p[j].rho,2)+epsilon))*gradKernel(p[i],p[j],0);
            aijy = -m*(p[i].p/(pow(p[i].rho,2)+epsilon) + p[j].p/(pow(p[j].rho,2)+epsilon))*gradKernel(p[i],p[j],1);
            p[i].ax += aijx;
            p[i].ay += aijy;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    } 
  }
  if(gamm!=0){
    //calculation of normal vector 
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
    p[i].vxh=p[i].vx+p[i].ax*dt/2;
    p[i].vyh=p[i].vy+p[i].ay*dt/2;
    p[i].vx+=p[i].ax*dt;
    p[i].vy+=p[i].ay*dt;
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
  }
}

void boundaryCondition(Particle_State p[])
{
  double LW=0.05*3;
  double RW=4.0-0.05*2;
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
    p[i].vx=p[i].vxh+p[i].ax*dt/2;
    p[i].vy=p[i].vyh+p[i].ay*dt/2;
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
  }
  
}


