#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"SPH.h"
#include"numbers.h"


void calcPsi(Particle_State p[], double Psi[], int bfst[], int nxt[])
{
  double delta[OBP];
  int i;

  for(i=0; i<OBP; i++){
    Psi[i]=0;
    delta[i]=0;
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
          int jb = 0;
          jb = jx + jy*nBx;
          int j = bfst[jb];
          //fprintf(stderr,"%d bfst accessed, %d %d %d\n", jb, i, jx, jy);
          if(j==-1){
	    continue;
	  }
          for(;;){
            if(j<FLP+BP){
              j=nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            delta[i-FLP-BP]+=kernel(p[i], p[j]);
            j=nxt[j];
            if(j==-1){
              break;
            }
          }
        }
      }
    }
  }

  for(i=0; i<OBP; i++){
    Psi[i]=rho0*delta[i];
  }
}

void AkinciCalcDensity(Particle_State p[], double Psi[], int bfst[], int nxt[])
{
  int i;
  for(i=0; i<N; i++){
    p[i].rho=0;
  }

  //  fprintf(stderr,"density initialized\n");
#pragma omp parallel for schedule(dynamic,64)
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
          //fprintf(stderr,"%d bfst accessed, %d %d %d\n", jb, i, jx, jy);
          if(j==-1){
	    continue;
	  }
          for(;;){
            double rhoij=0;
            if(j<FLP+BP){
            rhoij=p[i].mass*kernel(p[i], p[j]);
            }else if(i>=FLP+BP){
              rhoij=Psi[j-FLP-BP]*kernel(p[i], p[j]);
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
  //  fprintf(stderr, "Density calculated\n");

}

void AkinciCalcAccelByPressure(Particle_State p[], double Psi[], int bfst[], int nxt[])
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
          if(j==-1){
	    continue;
	  }
          for(;;){
            double aijx=0;
            double aijy=0;
            if(j<FLP+BP){
              aijx=-p[i].mass*p[j].mass*(p[j].p/(p[i].rho*p[j].rho+epsilon))*gradKernel(p[i],p[j],0);
              aijx=-p[i].mass*p[j].mass*(p[j].p/(p[i].rho*p[j].rho+epsilon))*gradKernel(p[i],p[j],0);
            }else if(j>=FLP+BP){
              aijx=-Psi[j-FLP-BP]*p[i].p*gradKernel(p[i],p[j],0)/(pow(p[i].rho,2)+epsilon);
              aijy=-Psi[j-FLP-BP]*p[i].p*gradKernel(p[i],p[j],1)/(pow(p[i].rho,2)+epsilon);
            }
            p[i].ax += aijx;
            p[i].ay += aijy;
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
            double dist = dx*dx+dy*dy;
            if(j<FLP+BP){
              viscCoef=2.0*nu*h*cs/(p[i].rho+p[j].rho);
              viscCoef=-viscCoef*(dot)/(dist*dist+0.01*h*h);
              if(time<DAMPTIME){
                viscCoef=viscCoef*damper;
              }
              if(dot<0){
                aijx = -p[j].mass*viscCoef*gradKernel(p[i], p[j], 0);
                aijy = -p[j].mass*viscCoef*gradKernel(p[i], p[j], 1);
                //fprintf(stderr, "dot=%f %f %f aijx=%f, aijy=%f\n",dot, viscCoef, gradKernel(p[i], p[j], 0),  aijx, aijy);
              }else if(dot>=0){
                aijx=0; 
                aijy=0;
              }
            }else if(j>=FLP+BP){//force from rigid body
              viscCoef=nu*h*cs/(2.0*p[i].rho+epsilon);
              viscCoef=-viscCoef*dot/(dist*dist+0.01*h*h);
              if(time<DAMPTIME){
                viscCoef=viscCoef*damper;
              }

              if(dot<0){
                aijx=-Psi[j-FLP-BP]*viscCoef*gradKernel(p[i],p[j],0);
                aijy=-Psi[j-FLP-BP]*viscCoef*gradKernel(p[i],p[j],1);
              }else if(dot>=0){
                aijx=0; 
                aijy=0;
              }
            }
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
}

void rigidBodyTimeIntegration(Particle_State p[], double *omega, FILE *fp, int time)
{
  double Fx, Fy;
  double cmx, cmy;
  double torque=0;
  double inertia=0;
  double anglarMomentum=0;
  int i;
  Fx=0;
  Fy=0;
  cmx=0, cmy=0;
  for(i=0; i<OBP; i++){
    Fx+=m*p[i+FLP+BP].ax;
    Fy+=m*p[i+FLP+BP].ay;
  }

  for(i=0; i<OBP; i++){//calculating center of mass
    cmx+=p[i].px/OBP;
    cmy+=p[i].py/OBP;
  }

  for(i=0; i<OBP; i++){//calculating torque
    double dx, dy;
    dx=0, dy=0;
    dx=p[i+FLP+BP].px-cmx;
    dy=p[i+FLP+BP].py-cmy;
    torque+=dx*Fy-dy*Fx;
  }

  for(i=0; i<OBP; i++){//calculating moment of inertia
    double dx=0, dy=0;
    double dist=0;
    dx=p[i+FLP+BP].px-cmx;
    dy=p[i+FLP+BP].py-cmy;
    dist=dx*dx+dy*dy;
    inertia+=rigidMass*dist;
  }

  anglarMomentum=inertia*(*omega);
  anglarMomentum+=(torque*dt);
  *omega=anglarMomentum/inertia;

  for(i=0; i<OBP; i++){
    double dx=0, dy=0;
    double ax, ay;
    dx=p[i+FLP+BP].px-cmx;
    dy=p[i+FLP+BP].py-cmy;
    ax=(Fx);
    ay=(Fy);
    p[i].vxh+=ax*dt-(*omega)*dy;
    p[i].vyh+=ay*dt+(*omega)*dx;
    p[i].vx=p[i].vxh+ax/2.0;
    p[i].vy=p[i].vyh+ay/2.0;
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
  }

  fprintf(fp, "%f %f %f %f %f %f %f \n", (double)(time*dt), cmx, cmy, 10.0, 1.0, *omega, inertia);
}