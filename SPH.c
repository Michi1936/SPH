#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"SPH.h"
#include"numbers.h"

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

double poly6(Particle_State p1, Particle_State p2)
{
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double coef=4.0/(M_PI*pow(h,8));
  double val=0;
  if(dist<=1){
    val=coef*pow(h*h-dist*dist,3);
  }else{
    val=0;
  }
  return val;
}

double gradSpikey(Particle_State p1, Particle_State p2, int axis)
{
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double coeff_x = -30.0*dx/(M_PI*pow(h,5)*dist+epsilon);
  double coeff_y = -30.0*dy/(M_PI*pow(h,5)*dist+epsilon);
  double val=0;
  if(axis==0){//x direction 
    if(dist<=h){
      val=coeff_x*pow(h-dist,2);
    }else{
      val=0;
    }
    if(axis==1){//x direction 
      if(dist<=h){
        val=coeff_y*pow(h-dist,2);
      }else{
        val=0;
      }
    }
  }
  return val;
  
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

double kernelVal(double q)
{
  double ans=0;
  if(0<=q && q<=1){
    ans=cubicSpline1(q);
  }else if(1<=q && q<=2){
    ans=cubicSpline2(q);
  }else if(q<0){

  }
  return ans;
}

double gradKernel(Particle_State p1, Particle_State p2, int axis)//calculate gradient of kernel
{
  
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double coeff_x = (dx)/(dist*h+epsilon);
  double coeff_y = (dy)/(dist*h+epsilon);
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

double forwardGradKernel(Particle_State p1, Particle_State p2, int axis)
{
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h;
  double coeff_x = (dx)/(dist*h+epsilon);
  double coeff_y = (dy)/(dist*h+epsilon);
  double dh = 0.01;
  double val=0;
  if(axis==0){//x direction 
    val=coeff_x*(kernelVal(q+(dh/2.0))-kernelVal(q-(dh/2.0)))/dh;
  }else if(axis==1){
    val=coeff_y*(kernelVal(q+(dh/2.0))-kernelVal(q-(dh/2.0)))/dh;
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

void calcDensity(Particle_State p[], int bfst[], int nxt[])
{
  int i;
  for(i=0; i<N; i++){
    p[i].rho=0;
  }

  //  fprintf(stderr,"density initialized\n");
#pragma omp parallel for schedule(dynamic,64)
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
            rhoij=p[j].mass*kernel(p[i], p[j]);

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


//Muller(2005) pressure model is used
void calcPressure(Particle_State p[])
{
  //set p=0
  int i;
  double coef;
  coef=(rho0*pow(cs,2))/7.0;
  for(i=0; i<N; i++){
    p[i].p=0;
  }
  
#pragma omp parallel for schedule(dynamic,64)
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

void calcAccelByExternalForces(Particle_State p[])
{
  int i;
#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){
    double aijx, aijy;
    aijx     = 0;
    aijy     = - g;//gravitational force
    p[i].ax += aijx;
    p[i].ay += aijy;
  }

#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP+BP; i<N; i++){
    double aijx, aijy;
    aijx     = 0;
    aijy     = - g;//gravitational force
    p[i].ax += aijx;
    p[i].ay += aijy;
  }
}

void calcAccelByPressure(Particle_State p[], int bfst[], int nxt[])
{
  int i;
#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){  
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
            aijx=-p[j].mass*((p[i].p/pow(p[i].rho,2.0)) + (p[j].p/pow(p[j].rho,2.0)))*gradKernel(p[i], p[j], 0);
            aijy=-p[j].mass*((p[i].p/pow(p[i].rho,2.0)) + (p[j].p/pow(p[j].rho,2.0)))*gradKernel(p[i], p[j], 1);
            p[i].ax += aijx;
            p[i].ay += aijy;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    } 
  }

  #pragma omp parallel for schedule(dynamic,64)
  for(i=FLP+BP; i<N; i++){  
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
            aijx=-p[j].mass*((p[i].p/pow(p[i].rho,2.0)) + (p[j].p/pow(p[j].rho,2.0)))*gradKernel(p[i], p[j], 0);
            aijy=-p[j].mass*((p[i].p/pow(p[i].rho,2.0)) + (p[j].p/pow(p[j].rho,2.0)))*gradKernel(p[i], p[j], 1);
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


void calcAccelByViscosity(Particle_State p[], int bfst[], int nxt[], int time)
    //Muller(2005) Weakly compressible SPH for free surface flow model is used.
{
  int i;

#pragma omp parallel for schedule(dynamic,64)
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
            aijx=0, aijy=0;
            double viscCoef=0;
            double dx = (p[i].px-p[j].px);
            double dy = (p[i].py-p[j].py);
            double dvx = (p[i].vx-p[j].vx);
            double dvy = (p[i].vy-p[j].vy);
            double dot = dx*dvx+dy*dvy;
            double dist = dx*dx+dy*dy;
            viscCoef=2.0*nu*h*cs/(p[i].rho+p[j].rho);
            viscCoef=-viscCoef*(dot)/(dist*dist+0.01*h*h);
	    if(time<DAMPTIME){
	      viscCoef=viscCoef*30.0;
	    }
            if(dot<0){
              aijx = -p[j].mass*viscCoef*gradKernel(p[i], p[j], 0);
              aijy = -p[j].mass*viscCoef*gradKernel(p[i], p[j], 1);
              //fprintf(stderr, "dot=%f %f %f aijx=%f, aijy=%f\n",dot, viscCoef, gradKernel(p[i], p[j], 0),  aijx, aijy);
            }else if(dot>=0){
              aijx=0; 
              aijy=0;
            }
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }
  
#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP+BP; i<N; i++){
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
            double dx = (p[i].px-p[j].px);
            double dy = (p[i].py-p[j].py);
            double dvx = (p[i].vx-p[j].vx);
            double dvy = (p[i].vy-p[j].vy);
            double dot = dx*dvx+dy*dvy;
            double dist = dx*dx+dy*dy;
            viscCoef=2.0*nu*h*cs/(p[i].rho+p[j].rho);
            viscCoef=-viscCoef*(dot)/(dist*dist+0.01*h*h);
	    if(time<DAMPTIME){
	      viscCoef=viscCoef*30.0;
	    }

            if(dot<0){
              aijx = -p[j].mass*viscCoef*gradKernel(p[i], p[j], 0);
              aijy = -p[j].mass*viscCoef*gradKernel(p[i], p[j], 1);
              //fprintf(stderr, "dot=%f %f %f aijx=%f, aijy=%f\n",dot, viscCoef, gradKernel(p[i], p[j], 0),  aijx, aijy);
            }else if(dot>=0){
              aijx=0; 
              aijy=0;
            }
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j = nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }

}


double surfaceTensionCoefficient(double r)//caluculating cohesion term 
{
  double coef=32.0/(M_PI*pow(h,9));
  if(2*r>h && r<=h){
    return (coef*pow(h-r,3.0)*pow(r,3.0));
  }
  else if(r>0 && 2*r<=h){
    return (coef*(2.0*pow(h-r,3.0)*pow(r,3.0) - pow(h,6.0)/64.0));
  }
  else {
    return 0;
  }
}

void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int nxt[])
{
  int i;
  
  double nx[FLP];
  double ny[FLP];

  for(i=0; i<FLP; i++){
    nx[i]=0;
    ny[i]=0;
  }

  for(i=0; i<FLP; i++){
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
            if(j>=FLP){
            j = nxt[j];
            if(j==-1){
              break;
            }
            continue;
            }
            nx[i]+=h*p[j].mass*gradKernel(p[i],p[j],0)/(p[j].rho+epsilon);
            ny[i]+=h*p[j].mass*gradKernel(p[i],p[j],1)/(p[j].rho+epsilon);
            
            j=nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }
  //calculation of cohesion force and surface area minimization
  //Akinci(2013) cohesion and surface tension model is used
  //  fprintf(stderr, "calculation of n is completed.\n");

#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){
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
            if(j>=FLP){
            j = nxt[j];
            if(j==-1)break;
            continue;
            }
            double aijx, aijy;
            double dx, dy, dist, Kij;
            aijx=0, aijy=0;
            dx = (p[i].px-p[j].px);
            dy = (p[i].py-p[j].py);
              
            dist = sqrt(dx*dx+dy*dy);
            Kij=2.0*rho0/(p[i].rho+p[j].rho+epsilon);

            aijx=-Kij*gamm*p[i].mass*(p[j].mass*surfaceTensionCoefficient(dist)*dx/(dist+epsilon) + (nx[i]-nx[j]));
            aijx=-Kij*gamm*p[i].mass*(p[j].mass*surfaceTensionCoefficient(dist)*dy/(dist+epsilon) + (ny[i]-ny[j]));
              
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
  
}

double boundaryGamma(Particle_State p1, Particle_State p2){
  double val=0;
  double dx = p1.px-p2.px;
  double dy = p1.py-p2.py;
  double dist = sqrt(dx*dx+dy*dy);
  double q =dist/h;

  if(0<q && q<(2.0/3.0)){
    val=2.0/3.0;
  }else if(2.0/3.0<q && q<1.0){
    val=(2.0*q-(2.0/3.0)*q*q);
  }else if(1.0<q && q<2.0){
    val=pow(2.0-q,2)/2.0;
  }else{
    val=0;
  }

  return 0.02*cs*cs*val/(dist+epsilon);
}

void calcAccelByBoundaryForce(Particle_State p[], int bfst[], int nxt[])//Boundary force in Monaghan(2005) model
{
  int i;

#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;

      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb=jx+jy*nBx;
          int j=bfst[jb];

          if(j==-1)continue;
          for(;;){
            //fprintf(stderr, "%d \n", j);
            if(j<FLP){
              j = nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            double aijx, aijy;
            double dx = p[i].px-p[j].px;
            double dy = p[i].py-p[j].py;
            double dist = sqrt(dx*dx+dy*dy);
            aijx=0, aijy=0;
	    aijx=(p[j].mass/(p[i].mass+p[j].mass))*boundaryGamma(p[i],p[j])*dx/(dist+epsilon);
	    aijy=(p[j].mass/(p[i].mass+p[j].mass))*boundaryGamma(p[i],p[j])*dy/(dist+epsilon);
	    //aijx=boundaryGamma(p[i], p[j])*dx/(2.0*dist+epsilon);
	    //aijy=boundaryGamma(p[i], p[j])*dy/(2.0*dist+epsilon);
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j=nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }


#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP+BP; i<N; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;

      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb=jx+jy*nBx;
          int j=bfst[jb];

          if(j==-1)continue;
          for(;;){
            //fprintf(stderr, "%d \n", j);
            if(j<FLP || FLP+BP<=j){
              j = nxt[j];
              if(j==-1){
                break;
              }
	      continue;
            }
            double aijx, aijy;
            double dx = p[i].px-p[j].px;
            double dy = p[i].py-p[j].py;
            double dist = sqrt(dx*dx+dy*dy);
            aijx=0, aijy=0;
	    aijx=(p[j].mass/(p[i].mass+p[j].mass))*boundaryGamma(p[i],p[j])*dx/(dist+epsilon);
	    aijy=(p[j].mass/(p[i].mass+p[j].mass))*boundaryGamma(p[i],p[j])*dy/(dist+epsilon);
		    // aijx=boundaryGamma(p[i], p[j])*dx/(2.0*dist+epsilon);
		    // aijy=boundaryGamma(p[i], p[j])*dy/(2.0*dist+epsilon);
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j=nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }

}

double adhesionCoefficient(Particle_State p1, Particle_State p2)
{
  double val=0;
  double dx = p1.px-p2.px;
  double dy = p1.py-p2.py;
  double dist = sqrt(dx*dx+dy*dy);

  if(2.0*dist>h && dist<=h){
    val=pow((-4.0*dist*dist/h)+6.0*dist-2.0*h, 1.0/4.0);
  }else{
    val=0;
  }

  return 0.007*val/pow(h,3.25);
}

void calcAccelByAdhesion(Particle_State p[], int bfst[], int nxt[])
{
  int i;
  double boundaryParticleDensity[BP+OBP];
  double psi[BP+OBP];
  double beta=-1.0;

  for(i=0; i<BP+OBP; i++){
    boundaryParticleDensity[i]=0;
    psi[i]=0;
  }
#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP; i<N; i++){//calculating number density of boundary particles 
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;

      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb=jx+jy*nBx;
          int j=bfst[jb];

          if(j==-1)continue;
          for(;;){
            //fprintf(stderr, "%d \n", j);
            if(j<FLP){
              j = nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            boundaryParticleDensity[i-FLP]+=kernel(p[i], p[j]);

            j=nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }

  for(i=0; i<BP+OBP; i++){
    psi[i]=rho0*boundaryParticleDensity[i];
  }

#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;

      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb=jx+jy*nBx;
          int j=bfst[jb];

          if(j==-1)continue;
          for(;;){
            //fprintf(stderr, "%d \n", j);
            if(j<FLP){
              j = nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            double aijx, aijy;
            double dx = p[i].px-p[j].px;
            double dy = p[i].py-p[j].py;
            double dist = sqrt(dx*dx+dy*dy);
            aijx=0, aijy=0;

            aijx=-beta*p[i].mass*psi[j-FLP]*adhesionCoefficient(p[i], p[j])*dx/(dist*epsilon);
            aijy=-beta*p[i].mass*psi[j-FLP]*adhesionCoefficient(p[i], p[j])*dy/(dist*epsilon);

            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j=nxt[j];
            if(j==-1)break;
          }
        }
      }
    }
  }
}

void rotateRigidBody(Particle_State p[], double angVel)
{
  double gx, gy;
  int i;
  gx=0, gy=0;


  for(i=FLP+BP; i<N; i++){//calculating center of mass
    gx+=p[i].px/OBP;
    gy+=p[i].py/OBP;
  }

  //#pragma omp parallel for schedule(dynamic,64)  
    for(i=FLP+BP; i<N; i++){
    double dx=p[i].px-gx;
    double dy=p[i].py-gy;

    p[i].vx+=(-angVel*dy);
    p[i].vy+=(angVel*dx);
  }

  fprintf(stderr,"%f Rigid body is rotated.\n",angVel);
}

/*
void rigidBodyCorrection(Particle_State p[], int bfst[], int nxt[], FILE *fp)
{
  int i;
  double rgx[2], rgy[2];
  double rgxDashed, rgyDashed;
  double rxDashed[OBP], ryDashed[OBP];
  double inertia;
  double theta;
  rgx[0]=0;
  rgy[0]=0;
  rgx[1]=0;
  rgy[1]=0;
  rgxDashed=0;
  rgyDashed=0;
  inertia=0;
  theta=0;
  for(i=0; i<OBP; i++){
    rxDashed[i]=0;
    ryDashed[i]=0;
  }
  
  for(i=FLP+BP; i<N; i++){//center of mass before time development
    rgx[0]+=p[i].prepx/(double)OBP;
    rgy[0]+=p[i].prepy/(double)OBP;
  }

  for(i=FLP+BP; i<N; i++){//center of mass after time development
    rgx[1]+=p[i].px/(double)OBP;
    rgy[1]+=p[i].py/(double)OBP;
  }

  for(i=FLP+BP; i<N; i++){//moment of inertia before time development
    double dx = fabs(p[i].px-rgx[0]);
    double dy = fabs(p[i].py-rgy[0]);
    double dist = sqrt(dx*dx+dy*dy);
    inertia+=p[i].mass*dist*dist;
  }

  for(i=FLP+BP; i<N; i++){//difference of position during time development
    rxDashed[i-FLP-BP]=p[i].px-p[i].prepx;
    ryDashed[i-FLP-BP]=p[i].py-p[i].prepy;
  }

  for(i=0; i<OBP; i++){
    rgxDashed+=rxDashed[i]/(double)OBP;
    rgyDashed+=ryDashed[i]/(double)OBP;
}
  
  for(i=FLP+BP; i<N; i++){//calculating rotation angle theta
    theta+=p[i].mass*( rxDashed[i-FLP-BP]*(p[i].prepy-rgy[1])-ryDashed[i-FLP-BP]*(p[i].prepx-rgx[1]))
           /(inertia+epsilon);
  }
  
  for(i=FLP+BP; i<N; i++){
    rxDashed[i-FLP-BP]=rgxDashed+(cos(theta)-1)*(p[i].prepx-rgx[1])+sin(theta)*(p[i].prepy-rgy[1]);
    ryDashed[i-FLP-BP]=rgyDashed+(-sin(theta))*(p[i].prepx-rgx[1])+(cos(theta)-1)*(p[i].prepy-rgy[1]);
  }

  for(i=FLP+BP; i<N; i++){
    p[i].vx=rxDashed[i-FLP-BP]/dt;
    p[i].vy=ryDashed[i-FLP-BP]/dt;
  }

  for(i=FLP+BP; i<N; i++){
    p[i].px=p[i].prepx+rxDashed[i-FLP-BP];
    p[i].py=p[i].prepy+ryDashed[i-FLP-BP];
  }
}
*/

void rigidBodyCorrection(Particle_State p[], FILE *fp, int time, double angVel){
  double gx, gy;
  double inertia;
  double qx[OBP], qy[OBP];
  double Tx, Ty;
  double Rot;
  int i;

  gx=0, gy=0;
  inertia=0;
  Tx=0, Ty=0;
  Rot=0;
  int threshold=0;

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
    inertia+=p[i].mass*(qx[i-FLP-BP]*qx[i-FLP-BP] + qy[i-FLP-BP]*qy[i-FLP-BP]);
  }

  for(i=FLP+BP; i<N; i++){//calculating translation velocity
    Tx+=p[i].vxh/OBP;
    Ty+=p[i].vyh/OBP;
  }


  if(time>threshold){
    for(i=FLP+BP; i<N; i++){//calculating anglar velocity
      Rot+=p[i].mass*(qx[i-FLP-BP]*p[i].vyh-qy[i-FLP-BP]*p[i].vxh)/(inertia+epsilon);
    }
  }else if(time<=threshold){
    Rot=angVel;
  }


  //  for(i=FLP+BP; i<N; i++){
    //     fprintf(stderr,"%d gx:%f gy:%f %f %f Tx:%f Ty:%f Rot:%f vx:%f vy:%f\n",i, gx, gy,qx[i-FLP-BP], qy[i-FLP-BP], Tx, Ty, Rot, p[i].vx, p[i].vy);
  //}

  for(i=FLP+BP; i<N; i++){
    p[i].vxh=Tx-Rot*qy[i-FLP-BP];
    p[i].vyh=Ty+Rot*qx[i-FLP-BP];
  }

  for(i=FLP+BP; i<N; i++){
    //    fprintf(stderr,"%d gx:%f gy:%f %f %f Tx:%f Ty:%f Rot:%f vx:%f vy:%f\n",i, gx, gy,qx[i-FLP-BP], qy[i-FLP-BP], Tx, Ty, Rot, p[i].vx, p[i].vy);
  }

  for(i=FLP+BP; i<N; i++){
     p[i].px=p[i].prepx+p[i].vxh*dt;
     p[i].py=p[i].prepy+p[i].vyh*dt;
  }

  gx=0;
  gy=0;

  for (i=FLP+BP; i<N; i++){//calculating center of mass
    gx+=p[i].px/OBP;
    gy+=p[i].py/OBP;
  }
  
  for (i=FLP+BP; i<N; i++){//calculating vector between center of mass and ith particle
    qx[i-FLP-BP]=p[i].px-gx;
    qy[i-FLP-BP]=p[i].py-gy;
  }

  for(i=FLP+BP; i<N; i++){//calculating inertia
    inertia+=p[i].mass*(qx[i-FLP-BP]*qx[i-FLP-BP] + qy[i-FLP-BP]*qy[i-FLP-BP]);
  }

  for(i=FLP+BP; i<N; i++){//calculating translation velocity
    Tx+=p[i].vx/OBP;
    Ty+=p[i].vy/OBP;
  }

  if(time>threshold){
    for(i=FLP+BP; i<N; i++){//calculating anglar velocity
      Rot+=p[i].mass*(qx[i-FLP-BP]*p[i].vy-qy[i-FLP-BP]*p[i].vx)/(inertia+epsilon);
    }
  }else if(time<=threshold){
    Rot=angVel;
  }

  for(i=FLP+BP; i<N; i++){
    p[i].vx=Tx-Rot*qy[i-FLP-BP];
    p[i].vy=Ty+Rot*qx[i-FLP-BP];
  }
  
  if(time%50==0){
  fprintf(fp, "%f %f %f %f %f %f %f\n", (double)(time*dt), gx, gy, Tx, Ty, Rot, inertia);
  }
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
  }

    for(i=FLP+BP; i<N; i++){
      p[i].vxh=p[i].vx+p[i].ax*dt/2.0;
      p[i].vyh=p[i].vy+p[i].ay*dt/2.0;
      p[i].vx+=p[i].ax*dt;
      p[i].vy+=p[i].ay*dt;
      p[i].prepx=p[i].px;
      p[i].prepy=p[i].py;
      p[i].px+=p[i].vxh*dt;
      p[i].py+=p[i].vyh*dt;
    }
}

void leapfrogStep(Particle_State p[], int time)
{
  int i;

  for(i=0; i<FLP; i++){
    p[i].vxh+=p[i].ax*dt;
    p[i].vyh+=p[i].ay*dt;
    p[i].vx=p[i].vxh+p[i].ax*dt/2.0;
    p[i].vy=p[i].vyh+p[i].ay*dt/2.0;
    if(p[i].inRegion==0){
      p[i].vxh =0; 
      p[i].vyh =0; 
      p[i].vx=0;
      p[i].vy=0;
    }
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
  }
  
  if(time>DAMPTIME){
    for(i=FLP+BP; i<N; i++){
      p[i].vxh+=p[i].ax*dt;
      p[i].vyh+=p[i].ay*dt;
      
      p[i].vx=p[i].vxh+p[i].ax*dt/2.0;
      p[i].vy=p[i].vyh+p[i].ay*dt/2.0;
      if(p[i].inRegion==0){
        p[i].vxh =0; 
        p[i].vyh =0; 
        p[i].vx=0;
        p[i].vy=0;
      }
      p[i].prepx=p[i].px;
      p[i].prepy=p[i].py;
      //    fprintf(stderr, "%f %f ", p[i].px, p[i].py);
      p[i].px+=p[i].vxh*dt;
      p[i].py+=p[i].vyh*dt;
      //  fprintf(stderr, "%f %f \n", p[i].px, p[i].py);
    }
  }
}

