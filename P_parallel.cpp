// Practice for parallel computing

#include <iostream>
#include <omp.h>
#include <cmath>
#include <time.h>
#include "P_parallel.h"

using namespace std;

const int N = 1E8;

int main(){
	clock_t t_st, t_end;
	double Time1, Time2;
	
	int *v_a, *v_b, *v_c;
	v_a = new int[N];
	v_b = new int[N];
	v_c = new int[N];
	
	
	for (int i=0; i<N; i++){
		v_a[i] = i;
		v_b[i] = i;	
	}
	
	// Call vector add1 function
	t_st = clock();
	Vec_add1(N, v_a, v_b, v_c);
	
	// sum all values in the vector c
	int sum=0;
	for (int i=0; i<N; i++){
		sum = sum + v_c[i];
	}
	cout << "sum1 = " << sum << endl;
	
	t_end = clock();
	
	Time1 = (t_end - t_st)/(double)CLOCKS_PER_SEC;
	cout << "Cal time for add1 (sec) = " << Time1 << endl;
		
	
	
	
	// Call vector add2 function
	t_st = clock();
	Vec_add2(N, v_a, v_b, v_c);
	
	// sum all values in the vector c
	sum=0;
	for (int i=0; i<N; i++){
		sum = sum + v_c[i];
	}
	cout << "sum2 = " << sum << endl;
	
	t_end = clock();
	
	Time1 = (t_end - t_st)/(double)CLOCKS_PER_SEC;
	cout << "Cal time for add1 (sec) = " << Time1 << endl;
	
	
	
	
	// Call vector add3 function
	t_st = clock();
	Vec_add3(N, v_a, v_b, v_c);
	
	// sum all values in the vector c
	sum=0;
	#pragma omp parallel for reduction(+:sum)
	for (int i=0; i<N; i++){
		sum = sum + v_c[i];
	}
	cout << "sum3 = " << sum << endl;
	
	t_end = clock();
	
	Time1 = (t_end - t_st)/(double)CLOCKS_PER_SEC;
	cout << "Cal time for add1 (sec) = " << Time1 << endl;
	
	
	// test
	int test = 5;
	#pragma omp parallel firstprivate(test)
	cout << "id: " << omp_get_thread_num() << ", test = " << test << endl;
	cout << "id: " << omp_get_thread_num() << ", test = " << test << endl;
	
	return 0;
}
