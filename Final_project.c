#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/// global variables
double Lx = 1.0;
double Ly = 1.0;
double Tb = 100.0;    // Temperature at bottom
double Tl = 50.0;     // Temperature at left
double Tr = 75.0;     // Temperature at right
double qt = 90.0;      // Heat transfer at top
double src = 0.0L;     // Source term in Poisson eq.
double err_cri = 1E-6;
double err_max;        // maximum error in the whole domain
int N_iter;

/// function definition
void Init(double *T, int Nx, int Ny, double dy);
void Point_Solver(double *T, int N_block);
void Cal_bc(double *T, int N_block);
void Save(double *T, int Nx, int Ny, double dx, double dy);

int main(){
    /// constant variables
    int N_block = 64;
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double *T;             // Temperature
    int i, j;

    T = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));


    Init(T, Nx, Ny, dy);        // Initialization
    Point_Solver(T, N_block);   // Iterative method (Point SOR) to calculate inner point
    Save(T, Nx, Ny, dx, dy);    // Save data


    return 0;
}



/// initialization
void Init(double *T, int Nx, int Ny, double dy){
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
    // Set top temperature
    for (i=1; i<Nx; i++)
        T[i+Ny*(Nx+1)] = T[i+(Ny-1)*(Nx+1)] - qt*dy;
}



/// Iterative method (Point SOR) to calculate inner point
void Point_Solver(double *T, int N_block){
    int Nx1 = N_block;
    int Nx2 = 2*N_block;
    int Nx = 3*N_block;
    int Ny1 = N_block;
    int Ny2 = 2*N_block;
    int Ny = 3*N_block;
    double dx = Lx/Nx;
    double dy = Ly/Ny;

    double *T_new;
    double Cp = -2.0L*(1.0L/dx/dx + 1.0L/dy/dy);
    double Ce = 1.0L/dx/dx;
    double Cw = 1.0L/dx/dx;
    double Cn = 1.0L/dy/dy;
    double Cs = 1.0L/dy/dy;
    double S = src;
    double err;
    int i, j;

    T_new = (double*) malloc((Nx+1)*(Ny+1)*sizeof(double));


    for (N_iter=0; N_iter<1E6; N_iter++){
        err_max = 0.0L;
        for (j=1; j<Ny; j++){
            // j outside the hole
            if (j<Ny1 || j>Ny2){
                for (i=1; i<Nx; i++){
                    T_new[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                                   -Cw*T[(i-1)+j*(Nx+1)] \
                                                   -Cn*T[i+(j+1)*(Nx+1)] \
                                                   -Cs*T[i+(j-1)*(Nx+1)] \
                                                   +S);
                    // Calculate err
                    err = fabs(T_new[i+j*(Nx+1)] - T[i+j*(Nx+1)]);
                    if (err>err_max)
                        err_max = err;
                    // Update new value
                    T[i+j*(Nx+1)] = T_new[i+j*(Nx+1)];
                }
            }
            // j inside the hole including BC
            else{
                for (i=1; i<Nx1; i++){
                    T_new[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                                   -Cw*T[(i-1)+j*(Nx+1)] \
                                                   -Cn*T[i+(j+1)*(Nx+1)] \
                                                   -Cs*T[i+(j-1)*(Nx+1)] \
                                                   +S);
                    // Calculate err
                    err = fabs(T_new[i+j*(Nx+1)] - T[i+j*(Nx+1)]);
                    if (err>err_max)
                        err_max = err;
                    // Update new value
                    T[i+j*(Nx+1)] = T_new[i+j*(Nx+1)];
                }
                for (i=Nx2+1; i<Nx; i++){
                    T_new[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                                   -Cw*T[(i-1)+j*(Nx+1)] \
                                                   -Cn*T[i+(j+1)*(Nx+1)] \
                                                   -Cs*T[i+(j-1)*(Nx+1)] \
                                                   +S);
                    // Calculate err
                    err = fabs(T_new[i+j*(Nx+1)] - T[i+j*(Nx+1)]);
                    if (err>err_max)
                        err_max = err;
                    // Update new value
                    T[i+j*(Nx+1)] = T_new[i+j*(Nx+1)];
                }
            }
        }
        /// Calculate new boundary values
        Cal_bc(T, N_block);

        /// Show iteration status
        printf("N_iter = %d, err_max = %g \n", N_iter, err_max);

        /// Consider if breaking the loop
        if (err_max<err_cri)
            break;
    }

    free(T_new);
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

    // Top outer bc
    for (i=1; i<Nx; i++)
        T[i+Ny*(Nx+1)] = T[i+(Ny-1)*(Nx+1)] - qt*dy;

    // bottom inner bc
    for (i=Nx1; i<=Nx2; i++){
        j = Ny1;
        if (i==Nx1){
            T[i+j*(Nx+1)] = 0.5L*( T[(i-1)+j*(Nx+1)] \
                                  +T[i+(j-1)*(Nx+1)]);
        }
        else if (i==Nx2){
            T[i+j*(Nx+1)] = 0.5L*( T[(i+1)+j*(Nx+1)] \
                                  +T[i+(j-1)*(Nx+1)]);
        }
        else{
            T[i+j*(Nx+1)] = T[i+(j-1)*(Nx+1)];
        }
    }
    // top inner bc
    for (i=Nx1; i<=Nx2; i++){
        j = Ny2;
        if (i==Nx1){
            T[i+j*(Nx+1)] = 0.5L*( T[(i-1)+j*(Nx+1)] \
                                  +T[i+(j+1)*(Nx+1)]);
        }
        else if (i==Nx2){
            T[i+j*(Nx+1)] = 0.5L*( T[(i+1)+j*(Nx+1)] \
                                  +T[i+(j+1)*(Nx+1)]);
        }
        else{
            T[i+j*(Nx+1)] = T[i+(j+1)*(Nx+1)];
        }
    }
    // left inner bc
    for (j=Ny1+1; j<Ny2; j++){
        i = Nx1;
        T[i+j*(Nx+1)] = T[(i-1)+j*(Nx+1)];
    }
    // right inner bc
    for (j=Ny1+1; j<Ny2; j++){
        i = Nx2;
        T[i+j*(Nx+1)] = T[(i+1)+j*(Nx+1)];
    }
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
