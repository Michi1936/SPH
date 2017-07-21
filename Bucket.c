#include"SPH.h"
#include"numbers.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void allocateBucket(int **bfst, int **blst, int **nxt)//Allocating memory for bfst, blst, nxt;
{
  *bfst=(int *)malloc(sizeof(int)*nBxy);
  *blst=(int *)malloc(sizeof(int)*nBxy);
  *nxt=(int *)malloc(sizeof(int)*N);
}

void checkParticle(Particle_State p[])
{
  int i;
  for(i=0; i<N; i++){
    if(p[i].px>MAX_X || p[i].px<MIN_X || p[i].py>MAX_Y || p[i].py<MIN_Y)
    {
      p[i].inRegion=0;//inRegion is 0 if particle is out of the calculation region
    }
    else{
      p[i].inRegion=1;
    }
  }
  //  fprintf(stderr, "particle checked\n");
}

void makeBucket(int *bfst, int *blst, int*nxt, Particle_State p[])//putting particles into backet
{
  int i;
  for(i=0; i<nBxy; i++){
    bfst[i]=-1;
  }
  for(i=0; i<nBxy; i++){
    blst[i]=-1;
  }
  for(i=0; i<N; i++){
    nxt[i]=-1;
  }
  //fprintf(stderr, "check bkt init\n");
  for(i=0; i<N; i++){
    if(p[i].inRegion==1){
      int ix = (int)((p[i].px - MIN_X)/BktLgth)+1;
      int iy = (int)((p[i].py - MIN_Y)/BktLgth)+1;
      int ib = ix+iy*nBx;//index of backet
      int j = blst[ib];
      blst[ib]=i;//particle i is last particle in backet ib
      if(j==-1){
        bfst[ib]=i;//if particle i is
      }else{
        nxt[j]=i;
      }
    }
  }
  //  fprintf(stderr, "Bucket is allocated\n");
}



void freeBucket(int *bfst, int *blst, int*nxt)
{
  free(bfst);
  free(blst);
  free(nxt);
}
