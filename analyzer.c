#include <stdio.h>
#include <stdlib.h>

#define FLP 40000

int main(int argc, char *argv[])
{
  int targetIndex;
  FILE *data;
  float posx[FLP];
  float posy[FLP];
  float vx[FLP];
  float vy[FLP];
  float vabs[FLP];
  float rho[FLP];
  float pressure[FLP];

  char imgName[128];
  char fName[256];
  char tmpStr[256];

  int date[2];
  int i;
  int spaceCount=0;
  int indexCount=0;
  int particleIndex=0;
  for(i=0; i<2; i++){
    date[i]=0;
  }

  
  fprintf(stderr, "Enter month-date:");
  scanf ("%d", &date[0]);
  fprintf(stderr, "Enter hour-minute:");
  scanf ("%d", &date[1]);
  fprintf(stderr, "Enter source image name:");
  scanf ("%s", imgName);

  fprintf(stderr, "Enter target index:\n");
  scanf("%d", &targetIndex);

  for(i=0; i<2; i++){
    fprintf(stderr, "%d\n", date[i]);
  }
  fprintf(stderr, "%s\n", imgName);

  sprintf(fName, "./Source_%s/%s_0.00_%d_%d_data.dat", imgName, imgName, date[0], date[1]);
  fprintf(stderr, "%s\n", fName);

  for(i=0; i<FLP; i++){
    posx[i]=0;
    posy[i]=0;
    vx[i]=0;
    vy[i]=0;
    vabs[i]=0;
    rho[i]=0; 
    pressure[i]=0;
  }
  fprintf(stderr, "data array initialized\n");

  data=fopen(fName, "r");
  if(data==NULL){
    fprintf(stderr, "File cannnot be opened!\n");
    exit(EXIT_FAILURE);
  }

  /* sprintf(fName, "./analyzerTest.dat"); */
  /* fprintf(stderr, "%s\n", fName); */
  /* data=fopen(fName, "r"); */

  if(data==NULL){
    fprintf(stderr, "File cannnot be opened!\n");
    exit(EXIT_FAILURE);
  }

  while(fgets(tmpStr, 256, data)!=NULL){
    if(*tmpStr=='\n'){
      spaceCount++;
      }

    if(indexCount==targetIndex && spaceCount==0){
      fprintf(stderr, "%s", tmpStr);
      sscanf(tmpStr, "%*d %f %f %f %f %f %f %f %*f %*f", &posx[particleIndex], &posy[particleIndex],
             &vx[particleIndex], &vy[particleIndex], &vabs[particleIndex],
             &rho[particleIndex], &pressure[particleIndex]);
    
      particleIndex++;
    }
 
     if(spaceCount==2){
      indexCount++;
      spaceCount=0;
     }
     
    if(indexCount>targetIndex){break;}

    //fprintf(stderr, "%s", tmpStr);

  }
  fprintf(stderr, "\nparticleIndex:%d index:%d\n", particleIndex,indexCount);
  for(i=0; i<3; i++){
    fprintf(stderr, "%f %f\n", posx[i], posy[i]);
}


  fclose(data);
  return 0;
}
