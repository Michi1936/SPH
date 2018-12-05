#include <stdio.h>
#include <stdlib.h>

#define FLP 40000

typedef struct{
  float px; 
  float py;
  float vx;
  float vy;
  float vabs;
  float ax;
  float ay;
  float rho; //density
  float pressure; //pressure
}pState;

void writePrevName(char imgName[], int date[]){
  FILE *fp;
  char *filename="prevDat.dat";
  char stringline[256];
  fp=fopen(filename, "w");
  sprintf(stringline, "%s %d %d", imgName, date[0], date[1]);
  fprintf(fp, "%s\n", stringline);
  fclose(fp);
}

void getPrevName(char imgName[], int date[]){
  FILE *fp;
  char *filename="prevDat.dat";
  char readline[256];
  fp=fopen(filename,"r");
  if(fp==NULL){
    fprintf(stderr, "prevDat.dat does not exist!\n");
    exit(0);
  }

  fgets(readline, 256, fp);
  sscanf(readline, "%s %d %d", imgName, &date[0], &date[1]);
  fprintf(stderr, "%s %d %d\n", imgName, date[0], date[1]);
  fclose(fp);
}

void inputFileNames(char imgName[], int date[]){
  fprintf(stderr, "Enter month-date:");
  scanf ("%d", &date[0]);
  fprintf(stderr, "Enter hour-minute:");
  scanf ("%d", &date[1]);
  fprintf(stderr, "Enter source image name:");
  scanf ("%s", imgName);
}

int main(void)
{
  int targetIndex;
  FILE *data;
  FILE *output;
  char imgName[128];
  char fName[256];
  char tmpStr[256];
  char outputfName[256];

  pState p[FLP];
  int date[2];
  int i;
  int spaceCount=0;
  int indexCount=0;
  int particleIndex=0;
  
  for(i=0; i<2; i++){
    date[i]=0;
  }
  for(i=0; i<FLP; i++){
    p[i].px=0;
    p[i].py=0;
    p[i].vx=0;
    p[i].vy=0;
    p[i].vabs=0;
    p[i].rho=0; 
    p[i].pressure=0;
  }

  fprintf(stderr, "1:enter file name\n");
  fprintf(stderr, "2:use previous settings\n");
  scanf("%d", &i);

  if(i==1){
    inputFileNames(imgName, date);
  }else if(i==2){
    getPrevName(imgName, date);
  }else{
    return 0;
  }

  fprintf(stderr, "Enter target index:\n");
  scanf("%d", &targetIndex);
 
 for(i=0; i<2; i++){
    fprintf(stderr, "%d\n", date[i]);
  }
  fprintf(stderr, "%s\n", imgName);

  //writing dat file name
  sprintf(fName, "./Source_%s/%s_0.00_%d_%d_data.dat", imgName, imgName, date[0], date[1]);
  fprintf(stderr, "%s\n", fName);

  writePrevName(imgName, date);

  //file opening
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

  //getting data at target index
  while(fgets(tmpStr, 256, data)!=NULL){
    if(*tmpStr=='\n'){
      spaceCount++;
    }

    //reach target index
    if(indexCount==targetIndex && spaceCount==0){
      fprintf(stderr, "%s", tmpStr);
      sscanf(tmpStr, "%*d %f %f %f %f %f %f %f %*f %*f", &p[particleIndex].px, &p[particleIndex].py,
             &p[particleIndex].vx, &p[particleIndex].vy, &p[particleIndex].vabs,
             &p[particleIndex].rho, &p[particleIndex].pressure);
      particleIndex++;
    }

    if(spaceCount==2){
      indexCount++;
      spaceCount=0;
    }
    if(indexCount>targetIndex){break;}
  }
  fprintf(stderr, "\nparticleIndex:%d index:%d\n", particleIndex,indexCount);
  for(i=0; i<3; i++){
    fprintf(stderr, "%f %f\n", p[i].px, p[i].py);
  }

  sprintf(outputfName, "INDEX%d_%s_%d_%d_.dat", targetIndex, imgName, date[0], date[1]);
  output=fopen(outputfName,"w");
  for(i=0; i<particleIndex; i++){
    fprintf(output,"%d %f %f %f %f %f\n", i, p[i].px, p[i].py, p[i].vx, p[i].vy, p[i].vabs);
  }

  fclose(data);
  fclose(output);
  return 0;
}
