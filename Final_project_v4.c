#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/// global variables
double Lx = 1.0;
double Ly = 1.0;
double Tb = 100.0;    // Temperature at bottom
double Tl = 50.0;     // Temperature at left
double Tr = 50.0;     // Temperature at right
double Ti = 75.0;     // Inner temperature
double qt = 90.0;     // Heat transfer at top
double res_cri = 1E-6;
double Max_Steps = 1E7;


/// function definition
void Init(double *T, double *Src, double *Res, int N_block);
void GaussSeidel_Iter(double *T, double *Src, int N_block);
void Multigrid_Iter(double *T, double *Src, int N_block);
void Cal_bc(double *T, int N_block);
void Save(double *T, int Nx, int Ny, double dx, double dy);
double Residual(double *Res, double *T, double *Src, int N_block);
void Set_Zero(double *array, int size);

int main(){
    /// constant variables
    int N_block = 32;
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double *T, *Src, *Res;
    double r, r1;
    int i, j;

    T = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));
    Src = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));
    Res = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));

    Init(T, Src, Res, N_block);        // Initialization

    r = Residual(Res, T, Src, N_block);
    int N_iter=0;
    while(r>res_cri && N_iter<Max_Steps){
        //Multigrid_Iter(T, Src, N_block);
        GaussSeidel_Iter(T, Src, N_block);
        r1 = Residual(Res, T, Src, N_block);
        printf("N_iter = %d, res = %g, ratio = %g\n", N_iter, r1, r1/r);
        r = r1;
        N_iter++;
        //system("pause");
    }
    Cal_bc(T, N_block);
    printf("N_block = %d, unknowns = %d\n", N_block, (Nx+1)*(Ny+1));
    Save(T, Nx, Ny, dx, dy);    // Save data

    return 0;
}



/// initialization
void Init(double *T, double *Src, double *Res, int N_block){
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    int i, j;

    // Set initial guess
    for (j=1; j<=Ny; j++)
        for (i=1; i<=Nx; i++)
            T[i+j*(Nx+1)] = 0.0L;
    // Set bottom temperature
    for (i=0; i<=Nx; i++)
        T[i+0*(Nx+1)] = Tb;
    // Set left temperature
    for (j=1; j<=Ny; j++)
        T[0+j*(Nx+1)] = Tl;
    // Set right temperature
    for (j=1; j<=Ny; j++)
        T[Nx+j*(Nx+1)] = Tr;
    // Set inner temperature
    for (j=Ny1; j<=Ny2; j++)
        T[Nx1+j*(Nx+1)] = Ti;
    for (j=Ny1; j<=Ny2; j++)
        T[Nx2+j*(Nx+1)] = Ti;
    for (i=Nx1; i<=Nx2; i++)
        T[i+Ny1*(Nx+1)] = Ti;
    for (i=Nx1; i<=Nx2; i++)
        T[i+Ny2*(Nx+1)] = Ti;

    // Set source
    for (j=0; j<=Ny; j++)
        for (i=0; i<=Nx; i++){
            Src[i+j*(Nx+1)] = 0.0;
            Res[i+j*(Nx+1)] = 0.0;
        }
    for (i=1; i<Nx; i++)
        Src[i+(Ny-1)*(Nx+1)] += qt/dy;
}



/// Iterative method (Point SOR) to calculate inner point
void GaussSeidel_Iter(double *T, double *Src, int N_block){
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = Lx/Nx;
    double dy = Ly/Ny;

    double Cp = -2.0*(1.0/dx/dx + 1.0/dy/dy);
    double Ce = 1.0/dx/dx;
    double Cw = 1.0/dx/dx;
    double Cn = 1.0/dy/dy;
    double Cs = 1.0/dy/dy;
    int i, j;

    for (j=1; j<Ny; j++){
        // j outside the hole
        if (j==Ny-1){
            for (i=1; i<Nx; i++){
                T[i+j*(Nx+1)] = (1.0/(Cp+Cn))*(-Ce*T[(i+1)+j*(Nx+1)] \
                                               -Cw*T[(i-1)+j*(Nx+1)] \
                                               -Cs*T[i+(j-1)*(Nx+1)] \
                                               +Src[i+j*(Nx+1)]);
            }
        }
        else if (j<Ny1 || j>Ny2){
            for (i=1; i<Nx; i++){
                T[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                           -Cw*T[(i-1)+j*(Nx+1)] \
                                           -Cn*T[i+(j+1)*(Nx+1)] \
                                           -Cs*T[i+(j-1)*(Nx+1)] \
                                           +Src[i+j*(Nx+1)]);
            }
        }
        else{
            for (i=1; i<Nx1; i++){
                T[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                           -Cw*T[(i-1)+j*(Nx+1)] \
                                           -Cn*T[i+(j+1)*(Nx+1)] \
                                           -Cs*T[i+(j-1)*(Nx+1)] \
                                           +Src[i+j*(Nx+1)]);
            }
            for (i=Nx2+1; i<Nx; i++){
                T[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                           -Cw*T[(i-1)+j*(Nx+1)] \
                                           -Cn*T[i+(j+1)*(Nx+1)] \
                                           -Cs*T[i+(j-1)*(Nx+1)] \
                                           +Src[i+j*(Nx+1)]);
            }
        }
    }
}



/// Multigrid Iteration
void Multigrid_Iter(double *T, double *Src, int N_block){
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;

    double *Res, *en, *rn;
    int i, j;


    if (N_block == 2){
        int N_iter = 0;
        int Max_Steps = 1000;
        double r;
        double res_cre = 1E-6;

        Res = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));
        r = 0.0;
        r = Residual(Res, T, Src, N_block);
        while(r>res_cri && N_iter<Max_Steps){
            GaussSeidel_Iter(T, Src, N_block);
            r = Residual(Res, T, Src, N_block);
            N_iter++;
        }
	}
    else{
        Res = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));
        en = (double*) malloc((Nx/2+1)*(Ny/2+1)*sizeof(double));
        rn = (double*) malloc((Nx/2+1)*(Ny/2+1)*sizeof(double));
        Set_Zero(Res, (Nx+1)*(Ny+1));
        Set_Zero(en, (Nx/2+1)*(Ny/2+1));
        Set_Zero(rn, (Nx/2+1)*(Ny/2+1));
        for (j=0; j<=Ny; j++)
            for (i=0; i<=Nx; i++)
                Res[i+j*(Nx+1)] = 0.0;
        // Smoother
        for (i=0; i<2; i++){
            GaussSeidel_Iter(T, Src, N_block);
        }
        Residual(Res, T, Src, N_block);

        // Projection: project the residue to next grid
        for (j=1; j<Ny/2; j++){
            for (i=1; i<Nx/2; i++){
                if (i<Nx1/2 || i>Nx2/2 || j<Ny1/2 || j>Ny2/2){
                    //rn[i+j*(Nx/2+1)] = Res[(i*2)+(j*2)*(Nx+1)];
                    // nine point
                    rn[i+j*(Nx/2+1)] = Res[((i*2)-1)+((j*2)-1)*(Nx+1)] + 2.0*Res[(i*2)+((j*2)-1)*(Nx+1)] + Res[((i*2)+1)+((j*2)-1)*(Nx+1)] \
                                     + 2.0*Res[((i*2)-1)+(j*2)*(Nx+1)] + 4.0*Res[(i*2)+(j*2)*(Nx+1)] + 2.0*Res[((i*2)+1)+(j*2)*(Nx+1)] \
                                     + Res[((i*2)-1)+((j*2)+1)*(Nx+1)] + 2.0*Res[(i*2)+((j*2)+1)*(Nx+1)] + Res[((i*2)+1)+((j*2)+1)*(Nx+1)];
                    rn[i+j*(Nx/2+1)] = 0.0625*rn[i+j*(Nx/2+1)];
                    //
                }
            }
        }
        // Calculate residue in next grid
        Multigrid_Iter(en, rn, N_block/2);
        // Interpolation
        for (j=0; j<=Ny; j++){
            for (i=0; i<=Nx; i++){
                if (i<=Nx1 || i>=Nx2 || j<=Ny1 || j>=Ny2){
                    /// We don't want to create new array for En
                    /// Here Res functions as En
                    if (i%2==0 && j%2==0){
                        Res[i+j*(Nx+1)] = 0.0;  // Empty the Res
                        Res[i+j*(Nx+1)] = en[(i/2)+(j/2)*(Nx/2+1)];
                    }
                    else if (i%2==1 && j%2==0){
                        Res[i+j*(Nx+1)] = 0.0;  // Empty the Res
                        Res[i+j*(Nx+1)] = 0.5*en[((i-1)/2)+(j/2)*(Nx/2+1)] \
                                         +0.5*en[((i+1)/2)+(j/2)*(Nx/2+1)];
                    }
                    else if (i%2==0 && j%2==1){
                        Res[i+j*(Nx+1)] = 0.0;  // Empty the Res
                        Res[i+j*(Nx+1)] = 0.5*en[(i/2)+((j-1)/2)*(Nx/2+1)] \
                                         +0.5*en[(i/2)+((j+1)/2)*(Nx/2+1)];
                    }
                    else if (i%2==1 && j%2==1){
                        Res[i+j*(Nx+1)] = 0.0;  // Empty the Res
                        Res[i+j*(Nx+1)] = 0.25*en[((i-1)/2)+((j-1)/2)*(Nx/2+1)] \
                                         +0.25*en[((i-1)/2)+((j+1)/2)*(Nx/2+1)] \
                                         +0.25*en[((i+1)/2)+((j-1)/2)*(Nx/2+1)] \
                                         +0.25*en[((i+1)/2)+((j+1)/2)*(Nx/2+1)];
                    }
                    else{
                        printf("Wrong!\n");
                    }
                }
            }
        }
        for (j=0; j<=Ny; j++){
            for (i=0; i<=Nx; i++){
                if (i<=Nx1 || i>=Nx2 || j<=Ny1 || j>=Ny2){
                    /// We don't want to create new array for En
                    /// Here Res functions as En
                    /// T adding the correction En
                    T[i+j*(Nx+1)] += Res[i+j*(Nx+1)];
                }
            }
        }
        for (i=0; i<2; i++){
            GaussSeidel_Iter(T, Src, N_block);
        }
        free(Res);
        free(en);
        free(rn);
    }
}



/// Calculate new boundary values
void Cal_bc(double *T, int N_block){
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    int i, j;

    // outer bc
    for (i=1; i<Nx; i++)
        T[i+Ny*(Nx+1)] = T[i+(Ny-1)*(Nx+1)] - qt*dy;
}



/// Save results
void Save(double *T, int Nx, int Ny, double dx, double dy){
    double x, y;
    int i, j;
    FILE *pFile, *pFile2;

    pFile2 = fopen("grid.txt","w");
    fprintf(pFile2, "%d\t%d", Nx, Ny);
    fclose(pFile2);

    pFile = fopen("Temperature.txt","w");
    for (j=0; j<=Ny; j++){
        for (i=0; i<=Nx; i++){
            x = dx*i;
            y = dy*j;
            fprintf(pFile, "%g\t%g\t%g\t\n", x, y, T[i+j*(Nx+1)]);
        }
    }
    fclose(pFile);
}



/// Gauss Elimination
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
    double dx = Lx/Nx;
    double dy = Ly/Ny;

    double Cp = -2.0L*(1.0L/dx/dx + 1.0L/dy/dy);
    double Ce = 1.0L/dx/dx;
    double Cw = 1.0L/dx/dx;
    double Cn = 1.0L/dy/dy;
    double Cs = 1.0L/dy/dy;
    double res_max = 0.0;
    int i, j;

    for (j=1; j<Ny; j++){
        // j outside the hole
        if (j==Ny-1){
            for (i=1; i<Nx; i++){
                Res[i+j*(Nx+1)] = Src[i+j*(Nx+1)] \
                                 -(Cp+Cn)*T[i+j*(Nx+1)] \
                                 -Ce*T[(i+1)+j*(Nx+1)] \
                                 -Cw*T[(i-1)+j*(Nx+1)] \
                                 -Cs*T[i+(j-1)*(Nx+1)];
                if (fabs(Res[i+j*(Nx+1)]) > res_max)
                    res_max = fabs(Res[i+j*(Nx+1)]);
            }
        }
        else if (j<Ny1 || j>Ny2){
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



/// Set zeros
void Set_Zero(double *array, int size){
    int i;
    for (i=0; i<size; i++)
        array[i] = 0.0;
}
