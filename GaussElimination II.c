#include <stdio.h>
#include <stdlib.h>

void GaussEli(double **A, double *b, double *x, int N);
double Residual(double *Res, double *T, double *Src, int N_block);

int main(){
    int N_block = 2;
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
    double *e, *Res, *Src;
    int i, j, k, m;

    A = (double **) malloc( (Nx+1)*(Ny+1)*sizeof(double*) );
    for (i=0; i<(Nx+1)*(Ny+1); i++)
        A[i] = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double));
    T = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );
    b = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );
    e = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );
    Res = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );
    Src = (double *) malloc( (Nx+1)*(Ny+1)*sizeof(double) );
    memset(Src, 0.0, (Nx+1)*(Ny+1));
    memset(Res, 0.0, (Nx+1)*(Ny+1));
    // Set T
    // Set initial guess
    for (j=1; j<=Ny; j++)
        for (i=1; i<=Nx; i++)
            T[i+j*(Nx+1)] = 0.0L;
    // Set bottom temperature
    for (i=0; i<=Nx; i++)
        T[i+0*(Nx+1)] = 100.0;
    // Set left temperature
    for (j=1; j<=Ny; j++)
        T[0+j*(Nx+1)] = 50.0;
    // Set right temperature
    for (j=1; j<=Ny; j++)
        T[Nx+j*(Nx+1)] = 75.0;
    // Set top temperature
    for (i=1; i<Nx; i++)
        T[i+Ny*(Nx+1)] = T[i+(Ny-1)*(Nx+1)] - 90.0*dy;

    for (j=1; j<=Ny; j++)
        for (i=1; i<=Nx; i++)
            //printf("i = %d, j = %d, x = %f\n", i, j, T[i+j*(Nx+1)]);

    // Set A, b
    for (j=0; j<(Nx+1)*(Ny+1); j++)
        for (i=0; i<(Nx+1)*(Ny+1); i++)
            A[j][i] = 0.0;

    for (j=0; j<=Ny; j++){
        for (i=0; i<=Nx; i++){
            if (j==0){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (i==0){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (i==Nx){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                b[i+j*(Nx+1)] = 0.0;
            }
            else if (j==Ny){
                A[i+j*(Nx+1)][i+j*(Nx+1)] = 1.0;
                A[i+j*(Nx+1)][i+(j-1)*(Nx+1)] = -1.0;
                b[i+j*(Nx+1)] = 0.0;
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

    b[1+1*(Nx+1)] = -5400.0;
    b[2+1*(Nx+1)] = -3600.0;
    b[3+1*(Nx+1)] = -3600.0;
    b[4+1*(Nx+1)] = -3600.0;
    b[5+1*(Nx+1)] = -6300.0;
    b[1+2*(Nx+1)] = -1800.0;
    b[5+2*(Nx+1)] = -2700.0;
    b[1+3*(Nx+1)] = -1800.0;
    b[5+3*(Nx+1)] = -2700.0;
    b[1+4*(Nx+1)] = -1800.0;
    b[5+4*(Nx+1)] = -2700.0;
    b[1+5*(Nx+1)] = -1260.0;
    b[2+5*(Nx+1)] = 540.0;
    b[3+5*(Nx+1)] = 540.0;
    b[4+5*(Nx+1)] = 540.0;
    b[5+5*(Nx+1)] = -2160.0;


    Residual(Res, T, Src, 2);
    for (j=0; j<=Ny; j++)
        for (i=0; i<=Nx; i++)
            printf("i = %d, j = %d, x = %f\n", i, j, Res[i+j*(Nx+1)]);
printf("============================\n");

    GaussEli(A, b, e, (Nx+1)*(Ny+1));

    for (j=0; j<=Ny; j++)
        for (i=0; i<=Nx; i++)
            T[i+j*(Nx+1)] += e[i+j*(Nx+1)];

    Residual(Res, T, Src, 2);


    for (j=0; j<=Ny; j++)
        for (i=0; i<=Nx; i++)
            printf("i = %d, j = %d, x = %f\n", i, j, Res[i+j*(Nx+1)]);

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


double Residual(double *Res, double *T, double *Src, int N_block){
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = 1.0/Nx;
    double dy = 1.0/Ny;

    double Cp = -2.0L*(1.0L/dx/dx + 1.0L/dy/dy);
    double Ce = 1.0L/dx/dx;
    double Cw = 1.0L/dx/dx;
    double Cn = 1.0L/dy/dy;
    double Cs = 1.0L/dy/dy;
    double res_max = 0.0;
    int i, j;

    for (j=1; j<Ny; j++){
        // j outside the hole
        if (j<Ny1 || j>Ny2){
            for (i=1; i<Nx; i++){
                Res[i+j*(Nx+1)] = Src[i+j*(Nx+1)] \
                                 -Cp*T[i+j*(Nx+1)] \
                                 -Ce*T[(i+1)+j*(Nx+1)] \
                                 -Cw*T[(i-1)+j*(Nx+1)] \
                                 -Cn*T[i+(j+1)*(Nx+1)] \
                                 -Cs*T[i+(j-1)*(Nx+1)];
                if (fabs(Res[i+j*(Nx+1)]) > res_max)
                    res_max = fabs(Res[i+j*(Nx+1)]);
            }
        }
        // j inside the hole including BC
        else{
            for (i=1; i<Nx1; i++){
                Res[i+j*(Nx+1)] = Src[i+j*(Nx+1)] \
                                 -Cp*T[i+j*(Nx+1)] \
                                 -Ce*T[(i+1)+j*(Nx+1)] \
                                 -Cw*T[(i-1)+j*(Nx+1)] \
                                 -Cn*T[i+(j+1)*(Nx+1)] \
                                 -Cs*T[i+(j-1)*(Nx+1)];
                if (fabs(Res[i+j*(Nx+1)]) > res_max)
                    res_max = fabs(Res[i+j*(Nx+1)]);
            }
            for (i=Nx2+1; i<Nx; i++){
                Res[i+j*(Nx+1)] = Src[i+j*(Nx+1)] \
                                 -Cp*T[i+j*(Nx+1)] \
                                 -Ce*T[(i+1)+j*(Nx+1)] \
                                 -Cw*T[(i-1)+j*(Nx+1)] \
                                 -Cn*T[i+(j+1)*(Nx+1)] \
                                 -Cs*T[i+(j-1)*(Nx+1)];
                if (fabs(Res[i+j*(Nx+1)]) > res_max)
                    res_max = fabs(Res[i+j*(Nx+1)]);
            }
        }
    }

    return res_max;
}
