#include <stdio.h>
#include <stdlib.h>

void GaussEli(double **A, double *b, double *x, int N);

int main(){
    int N_block = 8;
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = 1.0/Nx;
    double dy = 1.0/Ny;


    double **A;
    double *T;
    double *b;
    int i, j, k, m;

    A = (double **) malloc( (Nx+1)*(Ny+1)*sizeof(double*) );
    for (i=0; i<(Nx+1)*(Ny+1); i++)
        A[i] = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double));
    T = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );
    b = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );


    // Set A, b
    for (j=0; j<(Nx+1)*(Ny+1); j++)
        for (i=0; i<(Nx+1)*(Ny+1); i++)
            A[j][i] = 0.0;

    for (j=0; j<=Ny; j++){
        for (i=0; i<=Nx; i++){
            if (j==0){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 100.0;
            }
            else if (i==0){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 50.0;
            }
            else if (i==Nx){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 75.0;
            }
            else if (j==Ny){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][i+(j-1)*(Nx+1)] = -1.0;
                b[i+j*(Nx+1)] = -90.0*dy;
            }
            else if (j==Ny1 && i==Nx1){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][(i-1)+j*(Nx+1)] = -0.5;
                A[i+j*(Nx+1)][i+(j-1)*(Nx+1)] = -0.5;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (j==Ny1 && i==Nx2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][(i+1)+j*(Nx+1)] = -0.5;
                A[i+j*(Nx+1)][i+(j-1)*(Nx+1)] = -0.5;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (j==Ny2 && i==Nx1){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][(i-1)+j*(Nx+1)] = -0.5;
                A[i+j*(Nx+1)][i+(j+1)*(Nx+1)] = -0.5;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (j==Ny2 && i==Nx2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][(i+1)+j*(Nx+1)] = -0.5;
                A[i+j*(Nx+1)][i+(j+1)*(Nx+1)] = -0.5;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (j==Ny1 && i>Nx1 && i<Nx2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][i+(j-1)*(Nx+1)] = -1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (j==Ny2 && i>Nx1 && i<Nx2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][i+(j+1)*(Nx+1)] = -1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (i==Nx1 && j>Ny1 && j<Ny2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][(i-1)+j*(Nx+1)] = -1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (i==Nx2 && j>Ny1 && j<Ny2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][(i+1)+j*(Nx+1)] = -1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (i>Nx1 && i<Nx2 && j>Ny1 && j<Ny2){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else{
                A[i+j*(Nx+1)][i+j*(Nx+1)] = -2.0*(1.0/dx/dx + 1.0/dy/dy);
                A[i+j*(Nx+1)][(i-1)+j*(Nx+1)] = 1.0/dx/dx;
                A[i+j*(Nx+1)][(i+1)+j*(Nx+1)] = 1.0/dx/dx;
                A[i+j*(Nx+1)][i+(j-1)*(Nx+1)] = 1.0/dy/dy;
                A[i+j*(Nx+1)][i+(j+1)*(Nx+1)] = 1.0/dy/dy;
                b[i+j*(Nx+1)] = 0.0;
            }
        }
    }


    GaussEli(A, b, T, (Nx+1)*(Ny+1));

    for (i=0; i<(Nx+1)*(Ny+1); i++)
        printf("x = %f\n", T[i]);

    return 0;
}

void GaussEli(double **A, double *b, double *x, int N){
    int i, j, k, m;

    double temp = 0.0;
    /// Foreward elimination
    for(i=0; i<N-1; i++) {
        for(j=i+1; j<N; j++) {
            temp = A[j][i]/A[i][i];
            b[j] = b[j] - temp*b[i];
            for(k=i; k<N; k++) {
                A[j][k] = A[j][k] - temp*A[i][k];
            }
        }
    }
    /// Backward substitution
    double sum = 0.0;
    for(i=N-1; i>=0; i--) {
        sum = 0.0L;
        for(m=i+1; m<N; m++) {
            sum = sum + A[i][m]*x[m];
        }
        x[i] = (b[i]-sum)/A[i][i];
    }
}
