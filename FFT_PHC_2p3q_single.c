#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DEGUB 1

/// Function declaration
void bit_reverse(double *x_re, double *x_im, int N);
void butterfly3(double *x_re, double *x_im, int N_st, int N_end, int N);
void butterfly2(double *x_re, double *x_im, int N_st, int N_end, int N);

int main(){
    int i;
    const int N = 36;
    double *x_re, *x_im, *y_re, *y_im;
    int N_st, N_end;
    int d_2 = 0;
	int d_3 = 0; 
	int d_5 = 0;
    
    int n_temp = N;		// for calculation of each degree
    while(n_temp%2==0){
    	d_2++;
    	n_temp /= 2;
	}
    while(n_temp%3==0){
    	d_3++;
    	n_temp /= 3;
	}
    while(n_temp%5==0){
    	d_5++;
    	n_temp /= 5;
	}
    
    x_re = (double*) malloc(N*sizeof(double));
	x_im = (double*) malloc(N*sizeof(double));

    for (i=0; i<N; i++){
        x_re[i] = i;
        x_im[i] = 0;
    }

    // bit-reverse
    bit_reverse(x_re, x_im, N);
    
    // butterfly 3
    N_st = 1;
    N_end = pow(3,d_3);
    butterfly3(x_re, x_im, N_st, N_end, N);
    // butterfly 2
    N_st = N_end;
    N_end = N;
    butterfly2(x_re, x_im, N_st, N_end, N);
    
    // print results
    for (i=0; i<N; i++){
        printf("%f + %f i\n",x_re[i], x_im[i]);
    }
    
    return 0;
}


/// bit reverse for setting all input in place
void bit_reverse(double *x_re, double *x_im, int N){
	int base = 3;
    int m = N/base;    // added number which is highest bit
    int p, q;       // p & q is the index that exchanges for each other
    int k;          // k is use to check digital (log_2 k + 1) if is 1
    double temp;       // for temporary storage of number
    short *check;	// for check if no. has been changed 
    
    // set check array
    check = (short*) malloc(N*sizeof(short));
    int i;
	for (i=0; i<N; i++){
		check[i] = 0;
	}
	
	// for temp
	double *y_re, *y_im;
	y_re = (double*) malloc(N*sizeof(double));
	y_im = (double*) malloc(N*sizeof(double));
	for (i=0; i<N; i++){
		y_re[i] = x_re[i];
		y_im[i] = x_im[i];
	}
	
	
    // p is index before exchanging
    // from 0 to N-1 regularly
    // q is index after exchanging
    q = m;  // first index of exchanged number
    // skip p = 0 & N-1 because they are exchanged by themselves
    for (p=1; p<N-1; p++){
        #if DEGUB
        printf("%d <-> %d\n", p, q);
        #endif
        
        // for temp
        x_re[p] = y_re[q];
        x_im[p] = y_im[q];
                
        
        // exchange x to new place
        // if both no. have been changed, skip this exchange
        /*if (check[p]==0 || check[q]==0){
            temp = x_re[p];
            x_re[p] = x_re[q];
            x_re[q] = temp;
            temp = x_im[p];
            x_im[p] = x_im[q];
            x_im[q] = temp;
            
            // mark them as changed no.
            check[p] = 1;
            check[q] = 1;
        }*/

        // find next exchanged index q
        // add highest bit into old q
        k = m;
        base = 3;
        // if it needs to ¶i¦ì
        while(q>=(base-1)*k){
            q = q-(base-1)*k;    // (base-1) -> 0
            if (k%3 != 0)
            	base = 2;
            k = k/base;    // check next (right) digital
        }
        q = q+k;
    }
}


/// butterfly structure to recombine these input
void butterfly3(double *x_re, double *x_im, int N_st, int N_end, int N){
	int base = 3;
    int m;  // m is third of number in each group
            // it's also the difference between z and u
    int k;  // index in RHS of butterfly (output) (from 0 to m)
            // the rest of m/2 is counterpart of the first m/2
    int p;  // the first index in LHS of butterfly (input1)
    int q;  // the second index in LHS of butterfly (input2)
    int s;  // the third index in LHS of butterfly (input3)

    double w_k_re, w_k_im, w_N_re, w_N_im;
    double w_2k_re, w_2k_im, w_N_2_re, w_N_2_im;
    double w_but_re_3[3][3];
    double w_but_im_3[3][3];
    int i, j;
    for (i=1; i<base; i++){
    	for (j=1; j<base; j++){
    		w_but_re_3[i][j] =  cos((2.0D*M_PI/base)*i*j);
    		w_but_im_3[i][j] = -sin((2.0D*M_PI/base)*i*j);
		}
	}
	
    double temp, temp1_re, temp1_im, temp2_re, temp2_im, temp3_re;   // for temporary storage of number

    // loop for each steps of butterfly (step no.)
    for (m=N_st; m<N_end; m*=3){
        // Calculate multiplier of counterpart
        w_k_re = 1.0D;              // Re(W_{3^m}^0), multiplier of first counterpart
        w_k_im = 0.0D;              // Im(W_{3^m}^0)
        w_2k_re = 1.0D;             // Re(W_{3^2m}^0)
        w_2k_im = 0.0D;             // Re(W_{3^2m}^0)
        w_N_re =  cos(2.0D*M_PI/(m*3)); // Re(W_{3^m}^1)
        w_N_im = -sin(2.0D*M_PI/(m*3)); // Im(W_{3^m}^1)
        w_N_2_re = w_N_re*w_N_re - w_N_im*w_N_im;  // Re(W_{3^m}^2)
        w_N_2_im = 2.0D*w_N_re*w_N_im;  // Im(W_{3^m}^2)
                                // note that there is a minus sign for W_N = e^{-2PI/N}
        // loop for each output in a group (output no.)
        for (k=0; k<m; k++){
            // loop for each group (group no.)
            for (p=k; p<N; p+=3*m){
                // Calculate output of butterfly
                // find index of counterpart q
                q = p + m;
                s = q + m;
                // apply multiplier of counterpart
                // say, multiply W_{2^m}^k on x[q]
                temp = x_re[q];
                x_re[q] = w_k_re*x_re[q] - w_k_im*x_im[q];
                x_im[q] = w_k_re*x_im[q] + w_k_im*temp;
                temp = x_re[s];
                x_re[s] = w_2k_re*x_re[s] - w_2k_im*x_im[s];
                x_im[s] = w_2k_re*x_im[s] + w_2k_im*temp;

                // apply butterfly structure
                // to calculate x_p and x_q (counterpart).
                // here we calculate by input multiplied with
                // FFT_3_Matrix(multiplier on butterfly)
                temp1_re = x_re[p];
                temp2_re = x_re[q];
                temp3_re = x_re[s];
                x_re[p] = x_re[p] + x_re[q] + x_re[s];
                x_re[q] = temp1_re + (w_but_re_3[1][1]*x_re[q] - w_but_im_3[1][1]*x_im[q]) \
								   + (w_but_re_3[1][2]*x_re[s] - w_but_im_3[1][2]*x_im[s]);
                x_re[s] = temp1_re + (w_but_re_3[1][2]*temp2_re - w_but_im_3[1][2]*x_im[q]) \
								   + (w_but_re_3[2][2]*x_re[s] - w_but_im_3[2][2]*x_im[s]);
                temp1_im = x_im[p];
                temp2_im = x_im[q];
                x_im[p] = x_im[p] + x_im[q] + x_im[s];
                x_im[q] = temp1_im + (w_but_re_3[1][1]*x_im[q] + w_but_im_3[1][1]*temp2_re) \
								   + (w_but_re_3[1][2]*x_im[s] + w_but_im_3[1][2]*temp3_re);
                x_im[s] = temp1_im + (w_but_re_3[1][2]*temp2_im + w_but_im_3[1][2]*temp2_re) \
								   + (w_but_re_3[2][2]*x_im[s] + w_but_im_3[2][2]*temp3_re);
            }
            // calculate multiplier of next counterpart (with index k)
            temp = w_k_re;
            w_k_re = w_k_re*w_N_re - w_k_im*w_N_im;
            w_k_im = temp  *w_N_im + w_k_im*w_N_re;
            temp = w_2k_re;
            w_2k_re = w_2k_re*w_N_2_re - w_2k_im*w_N_2_im;
            w_2k_im = temp   *w_N_2_im + w_2k_im*w_N_2_re;
        }
    }
}


/// butterfly structure to recombine these input
void butterfly2(double *x_re, double *x_im, int N_st, int N_end, int N){
    int m;  // m is half of number in each group
            // it's also the difference between z and u
    int k;  // index in RHS of butterfly (output) (from 0 to m)
            // the rest of m/2 is counterpart of the first m/2
    int p;  // the first index in LHS of butterfly (input1)
    int q;  // the second index in LHS of butterfly (input2)

    double w_re, w_im, w_N_re, w_N_im;
    double temp;   // for temporary storage of number

    // loop for each steps of butterfly (step no.)
    for (m=N_st; m<N_end; m*=2){
        // Calculate multiplier of counterpart
        w_re = 1.0D;            // Re(W_{2^m}^0), multiplier of first counterpart
        w_im = 0.0D;            // Im(W_{2^m}^0)
        w_N_re = cos(M_PI/m);   // Re(W_{2^m}^1)
        w_N_im = -sin(M_PI/m);  // Im(W_{2^m}^1)
                                // note that there is a minus sign for W_N = e^{-2PI/N}
        // loop for each output in a group (output no.)
        for (k=0; k<m; k++){
            // loop for each group (group no.)
            for (p=k; p<N; p+=2*m){
                // Calculate output of butterfly
                // find index of counterpart q
                q = p + m;
                // apply multiplier of counterpart
                // say, multiply W_{2^m}^k on x[q]
                temp = x_re[q];
                x_re[q] = w_re*x_re[q] - w_im*x_im[q];
                x_im[q] = w_re*x_im[q] + w_im*temp;

                // apply butterfly structure
                // to calculate x_p and x_q (counterpart).
                // here we calculate by input multiplied with
                // FFT_2_Matrix(multiplier on butterfly)
                // ,which is [1 1; 1 -1]
                temp = x_re[p];
                x_re[p] = x_re[p] + x_re[q];
                x_re[q] = temp    - x_re[q];
                temp = x_im[p];
                x_im[p] = x_im[p] + x_im[q];
                x_im[q] = temp    - x_im[q];
            }
            // calculate multiplier of next counterpart (with index k)
            temp = w_re;
            w_re = w_re*w_N_re - w_im*w_N_im;
            w_im = temp*w_N_im + w_im*w_N_re;
        }
    }
}













