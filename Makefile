SPH:main.c SPH.c SPH_misc.c  setting.c Bucket.c SPH.h numbers.h Parameters.h
	gcc -o SPH -fopenmp main.c SPH.c SPH_misc.c setting.c Bucket.c -lm
