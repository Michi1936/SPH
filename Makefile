SPH:main.c SPH.c setting.c Bucket.c SPH.h numbers.h
	gcc -o SPH -fopenmp main.c SPH.c setting.c Bucket.c -lm
