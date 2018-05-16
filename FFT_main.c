#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "FFT_PHC.h"

int main(){
    int i;
    int N = pow(2,5)*pow(3,5)*pow(5,5) + 0.5D;
    double *x_re, *x_im;
	clock_t t1, t2;
    double time;


    // set raw data
    x_re = (double*) malloc(N*sizeof(double));
	x_im = (double*) malloc(N*sizeof(double));
    for (i=0; i<N; i++){
        x_re[i] = i;
        x_im[i] = 0;
    }


	t1 = clock();

    FFT(x_re, x_im, N);

	t2 = clock();
    time = (t2 - t1)/(double)CLOCKS_PER_SEC;

	// print results
    for (i=0; i<N; i++){
        //printf("%f + %f i\n",x_re[i], x_im[i]);
    }



    // show detail
    printf("N = %d\nTime (sec) = %f\n", N, time);

    return 0;
}
