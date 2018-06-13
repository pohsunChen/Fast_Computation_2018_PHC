#include <iostream>
#include <cmath>

/// constant variables
int N_block = pow(2,1);
int Nx1 = N_block;
int Nx2 = 2*N_block;
int Nx = 3*N_block;
int Ny1 = N_block;
int Ny2 = 2*N_block;
int Ny = 3*N_block;

long double Lx = 1.0L;
long double Ly = 1.0L;
long double dx = Lx/Nx;
long double dy = Ly/Ny;

long double Tb = 100.0L;    // Temperature at bottom
long double Tl = 50.0L;     // Temperature at left
long double Tr = 75.0L;     // Temperature at right
long double qt = 1.0L;      // Heat transfer at top
long double *T, *T_new;             // Temperature
long double src = 0.0L;     // Source term in Poisson eq.
long double err_cri = 1E-6;
long double err_max;        // maximum error in the whole domain
int N_iter;

/// function definition
void Alloc_mem();
void Init();
void Point_Solver();
void Cal_Boundary();



using namespace std;

int main(){

    Alloc_mem();    // Allocate memory
    cout << T[8] << endl;
    Init();         // Initialization
    cout << T[8] << endl;
    Point_Solver();
    cout << T[8] << endl;
    return 0;
}



/// Allocate memory
void Alloc_mem(){
    T = new long double((Nx+1)*(Ny+1));
    T_new = new long double((Nx+1)*(Ny+1));
}



/// initialization
void Init(){
    // Set initial guess
    for (int j=1; j<Ny; j++)
        for (int i=1; i<Nx; i++)
            T[i+j*(Nx+1)] = 0.0L;
    // Set bottom temperature
    for (int i=0; i<Nx; i++)
        T[i+0*(Nx+1)] = Tb;
    // Set left temperature
    for (int j=1; j<Ny; j++)
        T[0+j*(Nx+1)] = Tl;
    // Set right temperature
    for (int j=1; j<Ny; j++)
        T[Nx+j*(Nx+1)] = Tr;
    // Set top temperature
    for (int i=1; i<Nx; i++)
        T[i+Ny*(Nx+1)] = T[i+(Ny-1)*(Nx+1)] - qt*dy;
}


/// Iterative method (Point SOR) to calculate inner point
void Point_Solver(){
    long double Cp = -2.0L*(1.0L/dx/dx + 1.0L/dy/dy);
    long double Ce = 1.0L/dx/dx;
    long double Cw = 1.0L/dx/dx;
    long double Cn = 1.0L/dy/dy;
    long double Cs = 1.0L/dy/dy;
    long double S = src;
    long double err;

    err_max = 0.0L;
    for (N_iter=0; N_iter<1E6; N_iter++){
        for (int j=1; j<Ny; j++){
            // j outside the hole
            if (j<Ny1 || j>Ny2){
                for (int i=1; i<Nx; i++){
                    T_new[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                                   -Cw*T[(i-1)+j*(Nx+1)] \
                                                   -Cn*T[i+(j+1)*(Nx+1)] \
                                                   -Cs*T[i+(j-1)*(Nx+1)] \
                                                   +S);
                    // Calculate err
                    err = abs(T_new[i+j*(Nx+1)] - T[i+j*(Nx+1)]);
                    if (err>err_max)
                        err_max = err;
                    // Update new value
                    T[i+j*(Nx+1)] = T_new[i+j*(Nx+1)];
                }
            }
            // j inside the hole including BC
            else{
                for (int i=1; i<Nx1; i++){
                    T_new[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                                   -Cw*T[(i-1)+j*(Nx+1)] \
                                                   -Cn*T[i+(j+1)*(Nx+1)] \
                                                   -Cs*T[i+(j-1)*(Nx+1)] \
                                                   +S);
                }
                for (int i=Nx2+1; i<Nx; i++){
                    T_new[i+j*(Nx+1)] = (1.0L/Cp)*(-Ce*T[(i+1)+j*(Nx+1)] \
                                                   -Cw*T[(i-1)+j*(Nx+1)] \
                                                   -Cn*T[i+(j+1)*(Nx+1)] \
                                                   -Cs*T[i+(j-1)*(Nx+1)] \
                                                   +S);
                    // Calculate err
                    err = abs(T_new[i+j*(Nx+1)] - T[i+j*(Nx+1)]);
                    if (err>err_max)
                        err_max = err;
                    // Update new value
                    T[i+j*(Nx+1)] = T_new[i+j*(Nx+1)];
                }
            }
        }
        // Calculate new boundary values
        Cal_Boundary();

        //cout << "N_iter = " << N_iter << ", err_max = " << err_max << endl;
        // Consider if breaking the loop
        if (err_max<err_cri)
            break;
    }
}



void Cal_Boundary(){
    int i, j;
    // Top outer bc
    for (int i=1; i<Nx; i++)
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









