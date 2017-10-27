#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define h 0.1
int main()
{
  fprintf(stderr,"%f\n", (15.0/(pow(h,2)*14.0*M_PI)));
  return 0;
}
