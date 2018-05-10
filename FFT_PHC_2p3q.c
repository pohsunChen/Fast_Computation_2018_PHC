#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define DEGUB 1

/// Function declaration
void bit_reverse(double *x_re, double *x_im, int N, int base);
void butterfly(double *x_re, double *x_im, int N);

int main(){
    int i;
    const int N = 27;
    int base;
    double x_re[N], x_im[N];

    for (i=0; i<N; i++){
        x_re[i] = i;
        x_im[i] = 0.0D;
    }

    // bit-reverse
    base = 3;
    bit_reverse(x_re, x_im, N, base);
    // butterfly
    butterfly(x_re, x_im, N);

    // print results
    for (i=0; i<N; i++){
        printf("%f + %f i\n",x_re[i], x_im[i]);
    }

    return 0;
}


/// bit reverse for setting all input in place
void bit_reverse(double *x_re, double *x_im, int N, int base){
    int m = N/base;    // added number which is highest bit
    int p, q;       // p & q is the index that exchanges for each other
    int k;          // k is use to check digital (log_2 k + 1) if is 1
    double temp;       // for temporary storage of number


    // p is index before exchanging
    // from 0 to N-1 regularly
    // q is index after exchanging
    q = m;  // first index of exchanged number
    // skip p = 0 & N-1 because they are exchanged by themselves
    for (p=1; p<N-1; p++){
        #if DEGUB
        printf("%d <-> %d\n", p, q);
        #endif
        // do half of array
        if (p<q){
            temp = x_re[p];
            x_re[p] = x_re[q];
            x_re[q] = temp;
            temp = x_im[p];
            x_im[p] = x_im[q];
            x_im[q] = temp;
        }

        // find next exchanged index q
        // add highest bit into old q
        k = m;
        // if it needs to ¶i¦ì
        while(q>=(base-1)*k){
            q = q-(base-1)*k;    // (base-1) -> 0
            k = k/base;    // check next (right) digital
        }
        q = q+k;
    }
}


/// butterfly structure to recombine these input
void butterfly(double *x_re, double *x_im, int N){
    int m;  // m is half of number in each group
            // it's also the difference between z and u
    int k;  // index in RHS of butterfly (output) (from 0 to m)
            // the rest of m/2 is counterpart of the first m/2
    int p;  // the first index in LHS of butterfly (input1)
    int q;  // the second index in LHS of butterfly (input2)
    int s;  // the third index in LHS of butterfly (input3)

    double w_k_re, w_k_im, w_N_re, w_N_im;
    double w_2k_re, w_2k_im, w_N_2_re, w_N_2_im;
    double w_but_re_11 =  cos(2.0D*M_PI/3.0D);
    double w_but_im_11 = -sin(2.0D*M_PI/3.0D);
    double w_but_re_12 =  cos((2.0D*M_PI/3.0D)*2.0D);
    double w_but_im_12 = -sin((2.0D*M_PI/3.0D)*2.0D);
    double w_but_re_14 =  cos((2.0D*M_PI/3.0D)*4.0D);
    double w_but_im_14 = -sin((2.0D*M_PI/3.0D)*4.0D);
    double temp, temp1_re, temp1_im, temp2_re, temp2_im, temp3_re;   // for temporary storage of number

    // loop for each steps of butterfly (step no.)
    for (m=1; m<N; m*=3){
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
                x_re[q] = temp1_re + (w_but_re_11*x_re[q] - w_but_im_11*x_im[q]) + (w_but_re_12*x_re[s] - w_but_im_12*x_im[s]);
                x_re[s] = temp1_re + (w_but_re_12*temp2_re - w_but_im_12*x_im[q]) + (w_but_re_14*x_re[s] - w_but_im_14*x_im[s]);
                temp1_im = x_im[p];
                temp2_im = x_im[q];
                x_im[p] = x_im[p] + x_im[q] + x_im[s];
                x_im[q] = temp1_im + (w_but_re_11*x_im[q] + w_but_im_11*temp2_re) + (w_but_re_12*x_im[s] + w_but_im_12*temp3_re);
                x_im[s] = temp1_im + (w_but_re_12*temp2_im + w_but_im_12*temp2_re) + (w_but_re_14*x_im[s] + w_but_im_14*temp3_re);
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
















