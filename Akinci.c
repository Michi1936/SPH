#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"SPH.h"
#include"Parameters.h"
#include"numbers.h"

void setBoundaryType(int boundaryType[])
{
  int i;
  for(i=0; i<BP; i++){
    boundaryType[i]=0;
  }
  for(i=BP; i<BP+OBP; i++){
    boundaryType[i]=1;
  }
}

void calcPsi(Particle_State p[], double Psi[], int bfst[], int nxt[], int boundaryType[])
{
  double delta[BP+OBP];
  int i;
  for(i=0; i<BP+OBP; i++){
    Psi[i]=0;
    delta[i]=0;
  }
  
  //calculating Psi which is contribution of boundary particles 
#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP; i<N; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb = 0;
          jb = jx + jy*nBx;
          int j = bfst[jb];
          if(j==-1){
	    continue;
	  }
          for(;;){
            if(j<FLP){
              j=nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            if(boundaryType[i-FLP]!=boundaryType[j-FLP]){//Psi is added between the same particles
              j=nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            delta[i-FLP]+=kernel(p[i], p[j]);

            j=nxt[j];
            if(j==-1){
              break;
            }
          }
        }
      }
    }
  }

  for(i=0; i<BP+OBP; i++){
    Psi[i]=rho0/(delta[i]+epsilon);
  }
}

void AkinciCalcDensity(Particle_State p[], double Psi[], int bfst[], int nxt[])
{
  int i;
  for(i=0; i<N; i++){
    p[i].rho=0;
  }

#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){
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
          if(j==-1){
	    continue;
	  }
          for(;;){
            double rhoij=0;
            if(j<FLP){
              rhoij=p[j].mass*kernel(p[i], p[j]);
            }else if(i>=FLP){//interaction between boundary particles
              rhoij=Psi[j-FLP]*kernel(p[i], p[j]);
            }
            p[i].rho+=rhoij;
            j = nxt[j];
            if(j==-1){
	      break;
	    }
          }
        }
      }
    }
  }
}

void AkinciCalcAccelByPressure(Particle_State p[], double Psi[], int bfst[], int nxt[])
{
  int i;

#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP; i++){//calculating acceleration of fluid particles caused by pressure
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb = jx + jy*nBx;
          int j = bfst[jb];
          if(j==-1){
	    continue;
	  }
          for(;;){
            double aijx=0;
            double aijy=0;
            if(j<FLP){
              aijx=-p[i].mass*p[j].mass*(p[j].p/(p[i].rho*p[i].rho+epsilon))*gradKernel(p[i],p[j],0);
              aijy=-p[i].mass*p[j].mass*(p[j].p/(p[i].rho*p[i].rho+epsilon))*gradKernel(p[i],p[j],1);
            }else if(j>=FLP){//interaction betwen boundary particles
              aijx=-p[i].mass*Psi[j-FLP]*p[i].p*gradKernel(p[i],p[j],0)/(pow(p[i].rho,2)+epsilon);
              aijy=-p[i].mass*Psi[j-FLP]*p[i].p*gradKernel(p[i],p[j],1)/(pow(p[i].rho,2)+epsilon);
            }
            aijx=aijx/(p[i].mass+epsilon);
            aijy=aijy/(p[i].mass+epsilon);
            p[i].ax += aijx;
            p[i].ay += aijy;
            if(j>=FLP+BP){//calculation between fluid and rigid body
              p[j].ax+=-aijx;
	      p[j].ay+=-aijy;
            }
            j = nxt[j];
            if(j==-1){
	      break;
	    }
          }
        }
      }
    } 
  }
  
  //interaction between rigid body and wall particles
#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP+BP; i<N; i++){//acceleration of rigid body particles applied from boundary particles
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb = jx + jy*nBx;
          int j = bfst[jb];
          if(j==-1){
	    continue;
	  }
          for(;;){
            if(j<FLP || j>=FLP+BP){
              j=nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            double aijx=0;
            double aijy=0;
            aijx=-p[i].mass*Psi[j-FLP]*p[i].p*gradKernel(p[i],p[j],0)/(pow(p[i].rho,2)+epsilon);
            aijy=-p[i].mass*Psi[j-FLP]*p[i].p*gradKernel(p[i],p[j],1)/(pow(p[i].rho,2)+epsilon);
            p[i].ax += aijx/(p[i].mass+epsilon);
            p[i].ay += aijy/(p[i].mass+epsilon);
            j = nxt[j];
            if(j==-1){
	      break;
	    }
          }
        }
      }
    }
  }
}

void AkinciCalcAccelByViscosity(Particle_State p[], double Psi[], int bfst[], int nxt[], int time)
{
  int i;
  double damper=1.0;
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
          if(j==-1){
	    continue;
	  }
          for(;;){
            double aijx, aijy;
            aijx=0, aijy=0;
            double viscCoef=0;
            double dx = (p[i].px-p[j].px);
            double dy = (p[i].py-p[j].py);
            double dvx = (p[i].vx-p[j].vx);
            double dvy = (p[i].vy-p[j].vy);
            double dot = dx*dvx+dy*dvy;
            double dist = sqrt(dx*dx+dy*dy);
            if(j<FLP){//calculation between fluids particles
              viscCoef=2.0*nu*h*cs/(p[i].rho+p[j].rho);
	      viscCoef=-viscCoef*(dot)/(dist*dist+0.01*h*h);
	      aijx = -p[i].mass*p[j].mass*viscCoef*gradKernel(p[i], p[j], 0);
	      aijy = -p[i].mass*p[j].mass*viscCoef*gradKernel(p[i], p[j], 1);
            }else if(j>=FLP){//force from rigid body
              viscCoef=nu*h*cs/(2.0*p[i].rho+epsilon);
	      viscCoef=-viscCoef*dot/(dist*dist+0.01*h*h);
	      aijx=-p[i].mass*Psi[j-FLP]*viscCoef*gradKernel(p[i],p[j],0);
	      aijy=-p[i].mass*Psi[j-FLP]*viscCoef*gradKernel(p[i],p[j],1);
            }
            if(time<DAMPTIME){
              viscCoef=viscCoef*10.0;
            }
            aijx=aijx/(p[i].mass+epsilon);
            aijy=aijy/(p[i].mass+epsilon);
            p[i].ax+=aijx;
	    p[i].ay+=aijy;
            if(j>=FLP+BP){
              p[j].ax+=-aijx;
	      p[j].ay+=-aijy;
            }
            j = nxt[j];
            if(j==-1){
	      break;
	    }
          }
        }
      }
    }
  }
  damper=10.0;
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
          if(j==-1){
	    continue;
	  }
          for(;;){
            if(j<FLP || j>=FLP+BP){
              j=nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            double aijx, aijy;
            aijx=0, aijy=0;
            double viscCoef=0;
            double dx = (p[i].px-p[j].px);
            double dy = (p[i].py-p[j].py);
            double dvx = (p[i].vx-p[j].vx);
            double dvy = (p[i].vy-p[j].vy);
            double dot = dx*dvx+dy*dvy;
            double dist = sqrt(dx*dx+dy*dy);
            viscCoef=2.0*nu*h*cs/(p[i].rho+p[j].rho);
            viscCoef=-viscCoef*(dot)/(dist*dist+0.01*h*h);
            if(dot<0){
              aijx=-p[i].mass*Psi[j-FLP]*viscCoef*gradKernel(p[i],p[j],0);
              aijy=-p[i].mass*Psi[j-FLP]*viscCoef*gradKernel(p[i],p[j],1);
              //fprintf(stderr, "dot=%f %f %f aijx=%f, aijy=%f\n",dot, viscCoef, gradKernel(p[i], p[j], 0),  aijx, aijy);
            }else if(dot>=0){
              aijx=0; 
              aijy=0;
            }
            if(time<DAMPTIME){
              viscCoef=viscCoef*10.0;
            }

            p[i].ax+=aijx/(p[i].mass+epsilon);
	    p[i].ay+=aijy/(p[i].mass+epsilon);
            j = nxt[j];
            if(j==-1){
	      break;
	    }
          }
        }
      }
    }
  }
}

void rigidBodyTimeIntegration(Particle_State p[], RigidBodyValues rigV, FILE *fp, int time)
{
  double Fx=0, Fy=0;
  double transVx=0, transVy=0;
  double cmx=0, cmy=0;
  double Mass=0;
  double torque=0;
  double inertia=0;
  int i;

  if(time>MOTION_START_TIME){
    for(i=FLP+BP; i<N; i++){
      Mass+=p[i].mass;
    }
  
    for(i=FLP+BP; i<N; i++){//calculating F_total
      Fx+=rigidMass*p[i].ax;
      Fy+=rigidMass*p[i].ay;
    }
  
    for(i=FLP+BP; i<N; i++){//calculation of translation velocity
      transVx+=p[i].vx/OBP;
      transVy+=p[i].vy/OBP;
    }
    Fy+=-Mass*g;

    for(i=FLP+BP; i<N; i++){//calculating center of mass
      cmx+=p[i].px/OBP;
      cmy+=p[i].py/OBP;
    }

    for(i=FLP+BP; i<N; i++){//calculating torque
      double dx, dy;
      dx=0, dy=0;
      dx=p[i].px-cmx;
      dy=p[i].py-cmy;
      torque+=dx*Fy-dy*Fx;
    }
  
    for(i=FLP+BP; i<N; i++){//calculating moment of inertia
      double dx=0, dy=0;
      double dist=0;
      dx=p[i].px-cmx;
      dy=p[i].py-cmy;
      dist=dx*dx+dy*dy;
      inertia+=rigidMass*dist;
    }

    rigV.omega+=(torque/(inertia+epsilon))*dt;
    rigV.angle+=rigV.omega*dt;
    transVx+=Fx*dt/Mass;
    transVy+=Fy*dt/Mass;

    for(i=FLP+BP; i<N; i++){
      double dx=0, dy=0;
      dx=p[i].px-cmx;
      dy=p[i].py-cmy;
      p[i].vx=transVx-rigV.omega*dy;
      p[i].vy=transVy+rigV.omega*dx;
      p[i].px+=p[i].vx*dt;
      p[i].py+=p[i].vy*dt;
    }
  }
  fprintf(fp, "%f %f %f %f %f %f %f %f %f %f\n", (double)(time*dt), cmx, cmy, transVx, transVy, rigV.omega, inertia, torque, Fx, Fy);
}

void EulerCromerTimeIntegration(Particle_State p[])
{
  int i;
  for(i=0; i<FLP; i++){
    p[i].vx+=p[i].ax*dt;
    p[i].vy+=p[i].ay*dt;
    p[i].px+=p[i].vx*dt;
    p[i].py+=p[i].vy*dt;
  }
}
