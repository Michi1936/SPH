SPH:main_Akinci.c SPH.c SPH_misc.c setting.c Bucket.c SPH.h numbers.h Akinci.c
	gcc -o SPH -fopenmp main_Akinci.c SPH.c SPH_misc.c setting.c Bucket.c Akinci.c -lm
