#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define DEBUG 0

/// Function declaration
void bit_reverse(double *x_re, double *x_im, double *y_re, double*y_im, int N);
void butterfly(double *x_re, double *x_im, int N);
void set_base(int *base, int N);
int cal_id_b_max(int N);

int main(){
    int i;
    const int N = pow(2,2)*pow(3,0)*pow(5,2);
    double *x_re, *x_im, *y_re, *y_im;
    clock_t t1, t2;
    double time;

    x_re = (double*) malloc(N*sizeof(double));
	x_im = (double*) malloc(N*sizeof(double));
	y_re = (double*) malloc(N*sizeof(double));
	y_im = (double*) malloc(N*sizeof(double));

    for (i=0; i<N; i++){
        x_re[i] = i;
        x_im[i] = 0.0D;
    }

    t1 = clock();
    // bit-reverse (2 -> 3 -> 5)
    bit_reverse(x_re, x_im, y_re, y_im, N);
    // butterfly (5 -> 3 -> 2)
    butterfly(y_re, y_im, N);
    t2 = clock();
    time = (t2 - t1)/(double)CLOCKS_PER_SEC;


    // print results
    #if DEBUG
    for (i=0; i<N; i++){
        printf("%f + %f i\n",y_re[i], y_im[i]);
    }
    #endif

    // show detail
    printf("N = %d\nTime (sec) = %f", N, time);

    return 0;
}


/// bit reverse for setting all input in place
void bit_reverse(double *x_re, double *x_im, double *y_re, double*y_im, int N){
	int id_b_max = cal_id_b_max(N);
	int *base;
	int id_b;   // index for base
    int m;    // highest bit
    int p, q;       // p & q is the index that exchanges for each other
    int k;          // k is use to check digital (log_2 k + 1) if is 1

    // set base
	base = (int*) malloc((id_b_max+1)*sizeof(int));
	set_base(base, N);

    // p = 0 and N-1 are unchanged
	y_re[0] = x_re[0];
    y_im[0] = x_im[0];
	y_re[N-1] = x_re[N-1];
    y_im[N-1] = x_im[N-1];

    // p is index before exchanging
    // from 0 to N-1 regularly
    // q is index after exchanging
    m = N/base[id_b_max];    // added number which is highest bit
    q = m;  // first index of exchanged number
    // skip p = 0 & N-1 because they are exchanged by themselves
    for (p=1; p<N-1; p++){
        #if DEBUG
        printf("%d <-> %d\n", p, q);
        #endif

        // assign x to new place
        y_re[p] = x_re[q];
        y_im[p] = x_im[q];


        // find next exchanged index q
        // add highest bit into old q
        k = m;
        id_b = id_b_max;
        // if it needs to add 1 to right digital
        while(q>=(base[id_b]-1)*k){
            q = q-(base[id_b]-1)*k;    // (base-1) -> 0
            if ((k%base[id_b]) != 0){
            	id_b--;
			}
            k = k/base[id_b];    // check next (right) digital
        }
        q = q+k;
    }
}



/// butterfly structure to recombine these input
void butterfly(double *x_re, double *x_im, int N){
	int id_b_max = cal_id_b_max(N);   // calculate id_b_max
	int *base;
	int id_b;	// base id: 0 for base2, 1 for base3, and 2 for base5
    int m;  // m is number in each input group
            // it's also the difference between z and u
    int k;  // index in RHS of butterfly (output) (from 0 to m)
            // the rest of m/2 is counterpart of the first m/2
    int *id;  // the index in LHS of butterfly (input1-5)

	double *w_k_re, *w_k_im, *w_N_re, *w_N_im;  // multiplier of counterpart
    double ***w_but_re;     // multiplier on butterfly
    double ***w_but_im;
    double temp, *x_re_temp, *x_im_temp;   // for temporary storage of number
    int i, j, t;	// index for small loop


	// set base
	base = (int*) malloc((id_b_max+1)*sizeof(int));
	set_base(base, N);
	// allocate memory
	id = (int*) malloc(base[id_b_max]*sizeof(int));
    w_k_re = (double*) malloc(base[id_b_max]*sizeof(double));
    w_k_im = (double*) malloc(base[id_b_max]*sizeof(double));
    w_N_re = (double*) malloc(base[id_b_max]*sizeof(double));
    w_N_im = (double*) malloc(base[id_b_max]*sizeof(double));
    x_re_temp = (double*) malloc(base[id_b_max]*sizeof(double));
    x_im_temp = (double*) malloc(base[id_b_max]*sizeof(double));
    w_but_re = (double***) malloc((id_b_max+1)*sizeof(double**));
    w_but_im = (double***) malloc((id_b_max+1)*sizeof(double**));
    for(i=0; i<(id_b_max+1); i++){
        w_but_re[i] = (double**) malloc(base[id_b_max]*sizeof(double*));
        w_but_im[i] = (double**) malloc(base[id_b_max]*sizeof(double*));
        for (j=0; j<base[id_b_max]; j++){
            w_but_re[i][j] = (double*) malloc(base[id_b_max]*sizeof(double));
            w_but_im[i][j] = (double*) malloc(base[id_b_max]*sizeof(double));
        }
    }

    // initialization
    for (i=0; i<(id_b_max+1); i++){
    	for (j=0; j<base[id_b_max]; j++){
    		for (t=0; t<base[id_b_max]; t++){
    			w_but_re[i][j][t] = 0.0D;
    			w_but_im[i][j][t] = 0.0D;
			}
			w_k_re[j] = 0.0D;
			w_k_im[j] = 0.0D;
			w_N_re[j] = 0.0D;
			w_N_im[j] = 0.0D;
		}
	}

    // set the multiplier on butterfly
    for (id_b=0; id_b<(id_b_max+1); id_b++){
    	for (i=1; i<base[id_b_max]; i++){
	    	for (j=1; j<base[id_b_max]; j++){
	    		w_but_re[id_b][i][j] =  cos((2.0D*M_PI/base[id_b])*i*j);
	    		w_but_im[id_b][i][j] = -sin((2.0D*M_PI/base[id_b])*i*j);
			}
		}
	}


	m = 1;
	id_b = id_b_max;
    // loop for each steps of butterfly (step no.)
    while (m<N){
        // Calculate multiplier of counterpart
        for (i=1; i<base[id_b]; i++){
        	// i=0 is useless
        	w_k_re[i] = 1.0D;              // Re(W_{base^(i*m)}^0), multiplier counterpart
        	w_k_im[i] = 0.0D;              // Im(W_{base^(i*m)}^0)
		}
		// i=0 is useless
		w_N_re[1] =  cos(2.0D*M_PI/(m*base[id_b])); // Re(W_{5^m}^1)
        w_N_im[1] = -sin(2.0D*M_PI/(m*base[id_b])); // Im(W_{5^m}^1)
		for (i=2; i<base[id_b]; i++){
			w_N_re[i] = w_N_re[i-1]*w_N_re[1] - w_N_im[i-1]*w_N_im[1]; // Re(W_{5^(i*m)}^1)
        	w_N_im[i] = w_N_re[i-1]*w_N_im[1] + w_N_im[i-1]*w_N_re[1]; // Im(W_{5^(i*m)}^1)
		}

        // loop for each output (excluded counterpart) in a group (output no.)
        for (k=0; k<m; k++){
            // loop for each output group (no. of output group = (no of input group)/base)
            for (id[0]=k; id[0]<N; id[0]+=base[id_b]*m){
                // Calculate output of butterfly

				///// find index of counterpart p=id[0], q=id[1], s, t, l/////
                for (i=1; i<base[id_b]; i++)
                	id[i] = id[i-1] + m;


				///// apply multiplier of counterpart/////
                // say, multiply W_{2^m}^k on x[q]
                for (i=1; i<base[id_b]; i++){
                	temp = x_re[id[i]];
	                x_re[id[i]] = w_k_re[i]*x_re[id[i]] - w_k_im[i]*x_im[id[i]];
	                x_im[id[i]] = w_k_re[i]*x_im[id[i]] + w_k_im[i]*temp;
				}


                ///// apply butterfly structure/////
                // to calculate x_p and x_q (counterpart).
                // here we calculate by input multiplied with
                // FFT_3_Matrix(multiplier on butterfly)
                for (i=0; i<base[id_b]; i++){
                	x_re_temp[i] = x_re[id[i]];
                	x_im_temp[i] = x_im[id[i]];
				}

                /// real part
                for (j=1; j<base[id_b]; j++)
                	x_re[id[0]] += x_re[id[j]];
            	// loop for each counterpart
            	for (i=1; i<base[id_b]; i++){
            		x_re[id[i]] = x_re_temp[0];
            		// loop for each input
					for(j=1; j<base[id_b]; j++){
	            		x_re[id[i]] += (w_but_re[id_b][i][j]*x_re_temp[j] - w_but_im[id_b][i][j]*x_im_temp[j]);
					}
				}

				/// imaginary part
				for (j=1; j<base[id_b]; j++)
                	x_im[id[0]] += x_im[id[j]];
            	// loop for each counterpart
            	for (i=1; i<base[id_b]; i++){
            		x_im[id[i]] = x_im_temp[0];
            		// loop for each input
					for(j=1; j<base[id_b]; j++){
	            		x_im[id[i]] += (w_but_re[id_b][i][j]*x_im_temp[j] + w_but_im[id_b][i][j]*x_re_temp[j]);
					}
				}
            }

            // calculate multiplier of next counterpart (with index k)
            for (i=1; i<base[id_b]; i++){
            	temp = w_k_re[i];
	            w_k_re[i] = w_k_re[i]*w_N_re[i] - w_k_im[i]*w_N_im[i];
	            w_k_im[i] = temp*w_N_im[i] + w_k_im[i]*w_N_re[i];
			}
        }

        // check and change the base
        m *= base[id_b];
		if ((N/m)%base[id_b] != 0)
			id_b--;
    }
}


int cal_id_b_max(int N){
    int id_b_max = 0;
	if (N%2 == 0){
        id_b_max++;
	}
	if (N%3 == 0){
		id_b_max++;
	}
	if (N%5 == 0){
		id_b_max++;
	}
	if (N%7 == 0){
		id_b_max++;
	}
	return --id_b_max;
}

void set_base(int *base, int N){
    int i=0;
	if (N%2 == 0){
		base[i] = 2;
		i++;
	}
	if (N%3 == 0){
		base[i] = 3;
		i++;
	}
	if (N%5 == 0){
		base[i] = 5;
		i++;
	}
	if (N%7 == 0){
		base[i] = 7;
		i++;
	}
}
