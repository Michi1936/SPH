#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"SPH.h"
#include"Parameters.h"
#include"numbers.h"

//Monaghan(2005) cubic spline is used 
double cubicSpline1(double q)//cubic spline for 0<=q<=1
{
  return (15.0/(pow(h_smooth,2)*14.0*M_PI))*(pow(2.0-q,3)-4.0*pow(1.0-q,3));
}

double cubicSpline2(double q)//cubic spline for 1<=q<=2
{
  return (15.0/(pow(h_smooth,2)*14.0*M_PI))*(pow(2.0-q, 3));
}

double kernel(Particle_State p1, Particle_State p2)//p1 is central particle
{
  double dx = fabs(p1.px-p2.px);
  double dy = fabs(p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h_smooth;
  double val=0;
  if(0.0<=q && q<=1.0){
    val=cubicSpline1(q);
  }
  else if(1.0<=q && q<=2.0){
    val=cubicSpline2(q);
  }
  else {
    val=0;
  }
  return val;
}

double gradKernel(Particle_State p1, Particle_State p2, int axis)//calculate gradient of kernel
{
  
  double dx = (p1.px-p2.px);
  double dy = (p1.py-p2.py);
  double dist = sqrt(dx*dx+dy*dy);
  double q = dist/h_smooth;
  double coeff_x = (dx)/(dist*h_smooth+epsilon);
  double coeff_y = (dy)/(dist*h_smooth+epsilon);
  double val=0;
  if(axis==0){//x direction 
    if(0.0<=q && q<=1.0){
      val = (15.0/(pow(h_smooth,2)*14.0*M_PI))*coeff_x*(12.0*pow(1.0-q,2)-3.0*pow(2.0-q,2));
    }
    else if(1.0<=q && q<=2.0){
      val = -(15.0/(pow(h_smooth,2)*14.0*M_PI))*coeff_x*(3.0*pow(2.0-q,2));
    }
    else {
      val = 0;
    }
  }else if(axis==1){//y direction 
    if(0.0<=q && q<=1.0){
      val = (15.0/(pow(h_smooth,2)*14.0*M_PI))*coeff_y*(12.0*pow(1.0-q,2)-3.0*pow(2.0-q,2));
    }
    else if(1.0<=q && q<=2.0){
      val = -(15.0/(pow(h_smooth,2)*14.0*M_PI))*coeff_y*(3.0*pow(2.0-q,2));
    }
    else {
      val = 0;
    }
  }else{
    val = 0;
  }
  return val;
}

void calcDensity(Particle_State p[], int bfst[], int nxt[])
{
  int i;
  for(i=0; i<N; i++){
    p[i].rho=0;
  }

#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<N; i++){
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
            double rhoij=0;
            rhoij=p[j].mass*kernel(p[i], p[j]);
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

void integrateDensity(Particle_State p[], int bfst[], int nxt[])
{
  int i;
#pragma omp parallel for schedule(dynamic, 64)
  for(i=0; i<N; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb=jx+jy*nBx;
          int j=bfst[jb];
          if(j==-1){
	    continue;
	  }
          for(;;){
            double dvx = p[i].vx-p[j].vx;
            double dvy = p[i].vy-p[j].vy;    
            double dot = 0;
            dot=dvx*gradKernel(p[i],p[j],0)+dvy*gradKernel(p[i],p[j],1);
            p[i].rho+=p[j].mass*dot*dt;
            j=nxt[j];
            if(j==-1){
	      break;
	    }
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
  for(i=0; i<N; i++){
    p[i].p=0;
  }

  coef=(rho0*pow(cs,2))/7.0;  
#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<FLP+BP; i++){
    p[i].p = coef*(pow(p[i].rho/rho0,7)-1.0);//Tait equation
    if(p[i].p<0){
      p[i].p=0;
    }
  }

  //Pressure calculation of rigid body particles
  coef=(rigidMassMultiplier*rho0*pow(cs,2))/7.0;
#pragma omp parallel for schedule(dynamic,64)
  for(i=FLP+BP; i<N; i++){
    p[i].p = coef*(pow(p[i].rho/(rigidMassMultiplier*rho0),7)-1.0);//Tait equation
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
  for(i=0; i<N; i++){
    if(i>=FLP && i<FLP+BP){//gravitational force is not calculated for boundary particle
      continue;
    }
    double aijx, aijy;
    aijx     = 0;
    aijy     = -g;//gravitational force
    p[i].ax += aijx;
    p[i].ay += aijy;
  }
}

void calcAccelByPressure(Particle_State p[], int bfst[], int nxt[])
{
  int i;
#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<N; i++){  
    if(i>=FLP && i<FLP+BP){//if i is wall particle, skip this calculation
      continue;
    }
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
            double presCoef=0;
            presCoef=((p[i].p)/(pow(p[i].rho,2.0)+epsilon))+(p[j].p/(pow(p[j].rho,2.0)+epsilon));
            aijx=-p[j].mass*presCoef*gradKernel(p[i],p[j],0);
            aijy=-p[j].mass*presCoef*gradKernel(p[i],p[j],1);
            p[i].ax+=aijx;
            p[i].ay+=aijy;
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

void calcAccelByViscosity(Particle_State p[], int bfst[], int nxt[], int time)
    //Muller(2005) Weakly compressible SPH for free surface flow model is used.
{
  int i;
  double damper=10.0;
#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<N; i++){
    if(i>=FLP && i<FLP+BP){//wall particles are not accelerated hence calcAccel is skipped
      continue;
    }
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
            double aijx, aijy;
            double viscCoef=0;
            double alpha_visc=0;
            double dx = (p[i].px-p[j].px);
            double dy = (p[i].py-p[j].py);
            double dvx = (p[i].vx-p[j].vx);
            double dvy = (p[i].vy-p[j].vy);
            double dot = dx*dvx+dy*dvy;
            double dist = sqrt(dx*dx+dy*dy);
            aijx=0, aijy=0;      
            alpha_visc=2.0*nu*h_smooth*cs/(p[i].rho+p[j].rho+epsilon);
            viscCoef=-alpha_visc*dot/(dist*dist+epsilon);
            //double dot=dx*gradKernel(p[i],p[j],0)+dy*gradKernel(p[i],p[j],1);
            //aijx=4.0*nu*p[j].mass*dot*dvx/((p[i].rho+p[j].rho)*(dist*dist+epsilon));
	    //aijy=4.0*nu*p[j].mass*dot*dvy/((p[i].rho+p[j].rho)*(dist*dist+epsilon));
	    if(time<DAMPTIME){
	      viscCoef=viscCoef*damper;
	    }
	    aijx = -p[j].mass*viscCoef*gradKernel(p[i],p[j],0);
	    aijy = -p[j].mass*viscCoef*gradKernel(p[i],p[j],1);
            p[i].ax+=aijx;
            p[i].ay+=aijy;
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

double surfaceTensionCoefficient(double r)//caluculating cohesion term 
{
  double coef=32.0/(M_PI*pow(h_smooth,9));
  if(2*r>h_smooth && r<=h_smooth){
    return (coef*pow(h_smooth-r,3.0)*pow(r,3.0));
  }
  else if(r>0 && 2*r<=h_smooth){
    return (coef*(2.0*pow(h_smooth-r,3.0)*pow(r,3.0) - pow(h_smooth,6.0)/64.0));
  }
  else {
    return 0;
  }
}

//Based on Versatile Interactions at Interfaces for SPHbased Simulations(2016)
void calcInterfacialForce(Particle_State p[], int bfst[], int nxt[], FILE *fp)
{
  int i;
  double interfacialForce=0;
#pragma omp parallel for schedule(dynamic, 64)
  for(i=0; i<FLP; i++){
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
          double aijx, aijy;
          double interCoeff=0;
          double enl=ENLARGEMENT;//enlargement coefficient
          aijx=0, aijy=0;
          if(j<FLP){
            //interaction between fluid particles
            //some kind of surface tension
            interCoeff=FLUID_INTERACTION;
          }else{//interaction between fluid and boundary particles
            //1 means surface is hydrophily, 2 means surface is hydrophoby
            if(p[j].color==1){
              interCoeff=HPHILY_INTERACTION;
              //interCoeff=FLUID_INTERACTION*cos(CONTACT_ANGLE);
            }else if(p[j].color==2){
              interCoeff=HPHOBY_INTERACTION;
              //interCoeff=FLUID_INTERACTION*cos(CONTACT_ANGLE);
            }else if(p[j].color==3){
              interCoeff=FLUID_INTERACTION*(cos(CONTACT_ANGLE)+1.0)/2.0;//improved constant angle based upon Yang(2017)
            }else if(p[j].color==0){
              interCoeff=0;
            }
          }
	  double dx = (p[i].px-p[j].px);
	  double dy = (p[i].py-p[j].py);
          double dist = sqrt(dx*dx+dy*dy);
          
          if(dist<=enl*h_smooth){
            aijx=interCoeff*p[j].mass*cos((3.0*M_PI)*dist/(2.0*enl*h_smooth))*dx/(dist+epsilon);
	    aijy=interCoeff*p[j].mass*cos((3.0*M_PI)*dist/(2.0*enl*h_smooth))*dy/(dist+epsilon);
          }else{
            aijx=0;
            aijy=0;
          }

          p[i].ax+=aijx;
          p[i].ay+=aijy;
          interfacialForce+=sqrt(aijx*aijx+aijy*aijy);
          
          //if a fluid particle is forced by rigid body particle,
          //this rigid body particle must be forced by the fluid particle because of Newton's 3rd law.
          /* if(j>=FLP+BP){ */
          /*   p[j].ax+=-aijx; */
          /*   p[j].ay+=-aijy; */
          /* } */
          j = nxt[j];
          if(j==-1){
	    break;
	  }
        }
      }
    }
  }
 
#pragma omp parallel for schedule(dynamic, 64)
  for(i=FLP+BP; i<N; i++){
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
        //consider interaction between rigid and fluid
        for(;;){
          double aijx, aijy;
          double interCoeff=0;
          double enl=ENLARGEMENT;//enlargement coefficient
          aijx=0, aijy=0;
          //interaction between fluid and boundary particles
          //1 means surface is hydrophily, 2 means surface is hydrophoby
          if(p[i].color==1){
            interCoeff=HPHILY_INTERACTION;
            //interCoeff=FLUID_INTERACTION*cos(CONTACT_ANGLE);
          }else if(p[i].color==2){
            interCoeff=HPHOBY_INTERACTION;
            //interCoeff=FLUID_INTERACTION*cos(CONTACT_ANGLE);
          }else if(p[i].color==3){
            interCoeff=FLUID_INTERACTION*(cos(CONTACT_ANGLE)+1.0)/2.0;//improved constant angle based upon Yang(2017)
          }else if(p[i].color==0){
            interCoeff=0;
          }
	  double dx = (p[i].px-p[j].px);
	  double dy = (p[i].py-p[j].py);
          double dist = sqrt(dx*dx+dy*dy);
          
          if(dist<=enl*h_smooth){
            aijx=interCoeff*p[j].mass*cos((3.0*M_PI)*dist/(2.0*enl*h_smooth))*dx/(dist+epsilon);
	    aijy=interCoeff*p[j].mass*cos((3.0*M_PI)*dist/(2.0*enl*h_smooth))*dy/(dist+epsilon);
          }else{
            aijx=0;
            aijy=0;
          }
          if(j>=FLP){
            aijx=0;
            aijy=0;
          }

          p[i].ax+=aijx;
          p[i].ay+=aijy;
          interfacialForce+=sqrt(aijx*aijx+aijy*aijy);
          //if a fluid particle is forced by rigid body particle,
          //this rigid body particle must be forced by the fluid particle because of Newton's 3rd law.
          j = nxt[j];
          if(j==-1){
	    break;
	  }
        }
      }
    }
  }
   fprintf(fp, "%f\n", interfacialForce);
}

void calcAccelBySurfaceTension(Particle_State p[], int bfst[], int nxt[])
    //This surface tension model is based upon M. Becker & M.Teschner, Weakly compressible SPH for free surface flows
{
  int i;
#pragma omp parallel for schedule(dynamic,64)
  for(i=0; i<N; i++){
    if(i>=FLP && i<FLP+BP){
      continue;
    }
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
            if(j>=FLP){
              j=nxt[j];
              if(j==-1){
                break;
              }
              continue;
            }
            double dx=0;
            double dy=0;
            double surfCoef=0;
            dx=p[i].px-p[j].px;
            dy=p[i].py-p[j].py;

            surfCoef=-kappa*p[j].mass*kernel(p[i], p[j])/(p[i].mass);
            p[i].ax+=surfCoef*dx;
            p[i].ay+=surfCoef*dy;
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

double boundaryGamma(Particle_State p1, Particle_State p2){
  double val=0;
  double dx = p1.px-p2.px;
  double dy = p1.py-p2.py;
  double dist = sqrt(dx*dx+dy*dy);
  double q =dist/h_smooth;

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
  for(i=0; i<N; i++){
    if(i>=FLP && i<FLP+BP){
      continue;
    }
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px-MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py-MIN_Y)/BktLgth)+1;
      int jx, jy;
      for(jx=ix-1; jx<=ix+1; jx++){
        for(jy=iy-1; jy<=iy+1; jy++){
          int jb=jx+jy*nBx;
          int j=bfst[jb];
          if(j==-1){
	    continue;
	  }
          for(;;){
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
            p[i].ax+=aijx;
            p[i].ay+=aijy;
            j=nxt[j];
            if(j==-1){
	      break;
	    }
          }
        }
      }
    }
  }
}

void rotateRigidBody(Particle_State p[], RigidPreValue rig[], double angVel)
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
    double dx=rig[i-FLP-BP].prepx-gx;
    double dy=rig[i-FLP-BP].prepy-gy;

    p[i].vxh+=(-angVel*dy);
    p[i].vyh+=(angVel*dx);
  }

  for(i=FLP+BP; i<N; i++){
    double dx=p[i].px-gx;
    double dy=p[i].py-gy;

    p[i].vx+=(-angVel*dy);
    p[i].vy+=(angVel*dx);
  }

  fprintf(stderr,"\n%.2f Rigid body is rotated.\n",angVel);
}

void rigidBodyCorrection(Particle_State p[], RigidPreValue rig[], FILE *fp, int time, double com[]){
  double gx, gy;
  double inertia;
  double qx[OBP], qy[OBP];
  double Tx, Ty;
  double Rot;
  double radius;
  int i;

  gx=0, gy=0;
  inertia=0;
  Tx=0, Ty=0;
  Rot=0;
  radius=0;
  
  if(time>MOTION_START_TIME || time==1){
    for(i=0; i<OBP; i++){
      qx[i]=0;
      qy[i]=0;
    }
  
    for (i=0; i<OBP; i++){//calculating center of mass
      gx+=rig[i].prepx/OBP;
      gy+=rig[i].prepy/OBP;
    }

    for (i=0; i<OBP; i++){//calculating vector between center of mass and ith particle
      qx[i]=rig[i].prepx-gx;
      qy[i]=rig[i].prepy-gy;
    }

    for(i=0; i<OBP; i++){//calculating inertia
      double dist2=qx[i]*qx[i] + qy[i]*qy[i];
      inertia+=p[i+FLP+BP].mass*(dist2);
    }

    for(i=FLP+BP; i<N; i++){//calculating translation velocity
      Tx+=p[i].vxh/OBP;
      Ty+=p[i].vyh/OBP;
    }

    for(i=FLP+BP; i<N; i++){//calculating anglar velocity
      Rot+=p[i].mass*(qx[i-FLP-BP]*p[i].vyh-qy[i-FLP-BP]*p[i].vxh)/(inertia+epsilon);
    }

    for(i=FLP+BP; i<N; i++){//recalculating half step velocity
      p[i].vxh=Tx-Rot*qy[i-FLP-BP];
      p[i].vyh=Ty+Rot*qx[i-FLP-BP];
    }

    for(i=FLP+BP; i<N; i++){
      p[i].px=rig[i-FLP-BP].prepx+p[i].vxh*dt;
      p[i].py=rig[i-FLP-BP].prepy+p[i].vyh*dt;
    }
    //until here correction of position is complete
    //------------------------------------------------------
    //fprintf(stderr,"Check!\n");
    gx=0;
    gy=0;
    inertia=0;
    Rot=0;
    Tx=0;
    Ty=0;

    for (i=FLP+BP; i<N; i++){//calculating center of mass
      gx+=p[i].px/OBP;
      gy+=p[i].py/OBP;
    }
  
    for (i=FLP+BP; i<N; i++){//calculating vector between center of mass and ith particle
      qx[i-FLP-BP]=p[i].px-gx;
      qy[i-FLP-BP]=p[i].py-gy;
    }

    for(i=0; i<OBP; i++){//calculating inertia
      double dist2=qx[i]*qx[i] + qy[i]*qy[i];
      inertia+=p[i+FLP+BP].mass*(dist2);
    }

    for(i=FLP+BP; i<N; i++){//calculating translation velocity
      Tx+=p[i].vx/OBP;
      Ty+=p[i].vy/OBP;
    }

    for(i=FLP+BP; i<N; i++){//calculating anglar velocity
      Rot+=p[i].mass*(qx[i-FLP-BP]*p[i].vy-qy[i-FLP-BP]*p[i].vx)/(inertia+epsilon);
    }

    for(i=FLP+BP; i<N; i++){//recalculating velocity
      p[i].vx=Tx-Rot*qy[i-FLP-BP];
      p[i].vy=Ty+Rot*qx[i-FLP-BP];
    }
  }else if(time > 1 && time <=MOTION_START_TIME){//calculation while rigid body is fixed
    for(i=0; i<OBP; i++){
      qx[i]=0;
      qy[i]=0;
    }
  
    for (i=0; i<OBP; i++){//calculating center of mass
      gx+=rig[i].prepx/OBP;
      gy+=rig[i].prepy/OBP;
    }

    for (i=0; i<OBP; i++){//calculating vector between center of mass and ith particle
      qx[i]=rig[i].prepx-gx;
      qy[i]=rig[i].prepy-gy;
    }

    for(i=0; i<OBP; i++){//calculating inertia
      double dist=qx[i]*qx[i] + qy[i]*qy[i];
      inertia+=p[i+FLP+BP].mass*(dist);
    }
  }

  gx=0;
  gy=0;

  for (i=FLP+BP; i<N; i++){//calculating center of mass
    gx+=p[i].px/OBP;
    gy+=p[i].py/OBP;
  }

  for(i=0;i<OBP; i++)  {
    double dist=qx[i]*qx[i] + qy[i]*qy[i];
    if(dist>radius){
      radius=dist;
    }
  }
  //time, position of center of mass, translational velocity, angular velocity, inertia, radius of rigid body, spin parameter(ratio of anglar velocity at contact point between rigid body and fluid, to fall velocity)
  if(time%100==0){
    fprintf(fp, "%f %f %f %f %f %f %f %f %f\n", (double)(time*dt), gx, gy, Tx, Ty, Rot, inertia, radius, radius*Rot/(Ty+epsilon));
  }
  //outputting center of mass
  com[0]=gx;
  com[1]=gy;
}

void leapfrogStart(Particle_State p[], RigidPreValue rig[])
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
    rig[i-FLP-BP].prepx=p[i].px;
    rig[i-FLP-BP].prepy=p[i].py;
    p[i].px+=p[i].vxh*dt;
    p[i].py+=p[i].vyh*dt;
  } 
}

void leapfrogStep(Particle_State p[], RigidPreValue rig[], int time)
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
  
  if(time>MOTION_START_TIME){
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
      rig[i-FLP-BP].prepx=p[i].px;
      rig[i-FLP-BP].prepy=p[i].py;
      p[i].px+=p[i].vxh*dt;
      p[i].py+=p[i].vyh*dt;
    }
  }
}

