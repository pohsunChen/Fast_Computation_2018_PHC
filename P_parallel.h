// function for P_parallel.cpp

void Vec_add1(int N, int *v_a, int *v_b, int *v_c){
	for (int i=0; i<N; i++){
		v_c[i] = v_a[i] + v_b[i];
	}
	std::cout << "Method 1: \n";
}


void Vec_add2(int N, int *v_a, int *v_b, int *v_c){
	#pragma omp parallel for num_threads(4)
	for (int i=0; i<N; i++){
		v_c[i] = v_a[i] + v_b[i];
	}
	std::cout << "Method 2: \n";
}

void Vec_add3(int N, int *v_a, int *v_b, int *v_c){
	int N_thread = 4;
	int N_per_thread = N/N_thread;
	
	#pragma omp parallel num_threads(N_thread)
	{
		int id_th = omp_get_thread_num();
		
		for (int i=N_per_thread*id_th; i<(id_th+1)*N_per_thread; i++){
			v_c[i] = v_a[i] + v_b[i];
		}
	}
	std::cout << "Method 3: \n";
}


