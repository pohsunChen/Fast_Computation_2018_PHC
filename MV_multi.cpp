// matrix multiplication (parallel)

#include <iostream>
#include <iomanip>
#include <omp.h>
#include <time.h>
#include <cmath>
#include <stdlib.h>

using namespace std;

int main(){
	double *A, *x, *b1, *b2;
	int N = 10000;
	clock_t t1, t2;
	double time1, time2;
	int N_thread = 4;
	
	omp_set_num_threads(N_thread);
	
	// allocate memory
	//A = new double*[N];
	A = new double[N*N];
	x = new double[N];
	b1 = new double[N];
	b2 = new double[N];
	
	t1 = clock();
	// assign values
	srand(time(NULL));
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			A[i*N+j] = rand()%100*0.01;
		}
		x[i] = rand()%100*0.01;
	}
	
	// Multiply
	for (int i=0; i<N; i++){
		b1[i] = 0.0L;
		for (int j=0; j<N; j++){
			b1[i] += A[i*N+j]*x[j];
		}
	}
	t2 = clock();
	time1 = (t2 - t1)/(double) CLOCKS_PER_SEC;
	
	
	/*// Show results
	cout << "results: " << endl;
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			cout << setw(3) << A[i*N+j] << "\t";
		}
		cout << "|\t" << x[i] << "|\t" << b1[i] << endl;
	}
	*/
	
	
	
	//========================
	// Parallel     ==========
	//========================
		t1 = clock();
	// assign values
	int N_per_thread = N/4;
	int N_fin = N_per_thread*4;
	int N_rest = N%4;
	#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		srand(time(NULL)/(tid+1));
		for (int i=tid*N_per_thread; i<(tid+1)*N_per_thread; i++){
			for (int j=0; j<N; j++){
				A[i*N+j] = rand()%100*0.01;
			}
			x[i] = rand()%100*0.01;
		}
	}
	// Set rest of values
	if (N_rest!=0)
	#pragma omp parallel num_threads(N_rest)
	{
		int tid = omp_get_thread_num();
		srand(time(NULL)/(tid+1.5));
		for (int j=0; j<N; j++){
			A[(tid+N_fin)*N+j] = rand()%100*0.01;
		}
		x[(tid+N_fin)] = rand()%100*0.01;
	}
	
	// Check values
	/*for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			cout << A[i*N+j] << "\t";
		}
		cout << "\b\b| " << x[i] << endl;
	}
	*/
	
	int chunk = 10;
	# pragma omp parallel for schedule(dynamic,chunk)
	for (int i=0; i<N; i++){
		b2[i] = 0;
		for (int j=0; j<N; j++){
			b2[i] += A[i*N+j]*x[j];
		}
	}
	t2 = clock();
	time2 = (t2 - t1)/(double) CLOCKS_PER_SEC;
	
	
	// Show results
	cout << "time: (sec))\n";
	cout << "single\tparallel" << endl;
	cout << setw(5) << time1 << setw(12) << time2 << endl;
	cout << "\nSpeed up: " << time1/time2 << endl;
	
	delete A, x, b1, b2;
	return 0;
}
