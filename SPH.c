#include"SPH.h"
#include"numbers.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
//Equations below are based upon "Simulations of single bubbles rising through viscous liquids using SPH" by Szewc, Pozorski, Minier ,2013

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
  double coef=21.0/(16.0*M_PI*h*h*h);
  if(q<2.0){
    return coef*(pow(1.0-q/2.0,4)*(2.0*q+1.0));
  }else{
    return 0;
  }

}

double gradKernel(Particle_State p1, Particle_State p2, int dir)//calculate gradient of kernel
{
  
  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double coef=21.0/(16.0*M_PI*h*h*h);
  double coef_x=-coef*(5.0*q/(dist*h+epsilon))*(p1.px-p2.px);
  double coef_y=-coef*(5.0*q/(dist*h+epsilon))*(p1.py-p2.py);
  
  if(dir==0){//x direction 
    if(q<2.0){
      return coef_x*pow(1.0-q/2.0,3);
    }else{
      return 0;
    }
  }
  else if(dir==1){//y direction 
    if(q<2.0){
      return coef_y*pow(1.0-q/2.0,3);
    }else{
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


void calcPressure(Particle_State p[])
{
  //set p=0
  int i;
  for(i=0; i<N; i++){
    p[i].p=0;
  }
  
  for(i=0; i<N; i++){
    p[i].p = k1*(pow(p[i].rho/rho0,7)-1.0);//Tait equation
    
  }
}

void calcAcceleration(Particle_State p[], int bfst[], int blst[], int nxt[])
{
  int i,j;

  //inverse of the particle volume(Hu and Adams,2006)
  double Theta[FLP+BP];
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

  for(i=0; i<5; i++){
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
          //          fprintf(stderr, "i=%i, %f %f Theta_i=%f, Theta_j=%f\n", i, p[i].mu, p[j].mu, Theta[i], Theta[j]);
          //          fprintf(stderr, "rx=%f, ry=%f, dist=%f, gradKernel_x=%f, gradKernel_y=%f\n",
          //                  rx, ry, dist, gradKernel(p[i],p[j],0), gradKernel(p[i],p[j],1));
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

      
//calculation about pressure   
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
          while(j!=-1){
            double aijx, aijy;
            double presCoef=0;
            aijx=0, aijy=0;
            presCoef=(p[i].p/(pow(Theta[i],2)+epsilon)+p[j].p/(pow(Theta[j],2)+epsilon))/(p[i].mass+epsilon);

            aijx=presCoef*gradKernel(p[i], p[j], 0);
            aijy=presCoef*gradKernel(p[i], p[j], 1);
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j = nxt[j];
          }
        }
      }
    }
  }
  
  /*
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
            aijx = -p[i].mass*(p[i].p/(pow(p[i].rho,2)+epsilon) + p[j].p/(pow(p[j].rho,2)+epsilon))*gradKernel(p[i],p[j],0);
            aijy = -p[i].mass*(p[i].p/(pow(p[i].rho,2)+epsilon) + p[j].p/(pow(p[j].rho,2)+epsilon))*gradKernel(p[i],p[j],1);
            p[i].ax += aijx;
            p[i].ay += aijy;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    } 
    }*/

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


