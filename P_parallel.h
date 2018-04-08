// function for P_parallel.cpp

void Vec_add1(int N, int *v_a, int *v_b, int *v_c){
	for (int i=0; i<N; i++){
		v_c[i] = v_a[i] + v_b[i];
	}
}


void Vec_add2(int N, int *v_a, int *v_b, int *v_c){
	#pragma omp parallel for num_threads(16)
	for (int i=0; i<N; i++){
		v_c[i] = v_a[i] + v_b[i];
	}
}



