#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define DEBUG 0

/// Function declaration
void bit_reverse(double *x_re, double *x_im, int N);
void butterfly(double *x_re, double *x_im, int N);

int main(){
    int i;
    const int N = pow(2,16);
    double x_re[N], x_im[N];
    clock_t t1, t2;
    double time;

    for (i=0; i<N; i++){
        x_re[i] = i;
        x_im[i] = 0.0D;
    }

    t1 = clock();
    
    // bit-reverse
    bit_reverse(x_re, x_im, N);
    // butterfly
    butterfly(x_re, x_im, N);
    
	t2 = clock();
    time = (t2 - t1)/(double)CLOCKS_PER_SEC;


    // print results
    #if DEBUG
    for (i=0; i<N; i++){
        printf("%f + %f i\n",x_re[i], x_im[i]);
    }
    #endif

    // show detail
    printf("N = %d\nTime (sec) = %f", N, time);

    return 0;
}


/// bit reverse for setting all input in place
void bit_reverse(double *x_re, double *x_im, int N){
    int m = N/2;    // added number which is highest bit
    int p, q;       // p & q is the index that exchanges for each other
    int k;          // k is use to check digital (log_2 k + 1) if is 1
    double temp;       // for temporary storage of number


    // p is index before exchanging
    // from 0 to N-1 regularly
    // q is index after exchanging
    q = m;  // first index of exchanged number
    // skip p = 0 & N-1 because they are exchanged by themselves
    for (p=1; p<N-1; p++){
        #if DEBUG
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
        // take N = 16 for example, k = m = N/2 = 8
        // (1110) = 14 >= (k=8), log_2 k + 1 = 4
        // (0110) = 6  >= (k=4), log_2 k + 1 = 3
        // (0010) = 2  >= (k=2), log_2 k + 1 = 2
        // (0000) = 0   < (k=1) ---> break
        // (add 1 in digital log(k)+1 = 1)binary ---> (0001) = 1

        // take N = 6 for example
        // (0110) = 6   < (8=k) ---> break
        // (add 1 in digital log(k)+1 = 4)binary ---> (1110) = 14
        k = m;
        while(q>=k){
            q = q-k;    // 1 -> 0
            k = k/2;    // check next (right) digital
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

    double w_re, w_im, w_N_re, w_N_im;
    double temp;   // for temporary storage of number

    // loop for each steps of butterfly (step no.)
    for (m=1; m<N; m*=2){
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


















