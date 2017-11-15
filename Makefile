SPH:main.c SPH.c SPH_misc.c  setting.c Bucket.c SPH.h numbers.h
	gcc -o SPH -fopenmp main.c SPH.c SPH_misc.c setting.c Bucket.c -lm
