# include <iostream>
# include <omp.h>

using namespace std;

int main(){
	int tid;
	
	#pragma omp parallel for private(tid)
	for (int i=0; i<4; i++){
		tid = omp_get_thread_num();
		cout << "i: " << i << "\t tid: " << tid << endl;
	}
	cout << "+++++++++++++++++++++++++" << endl;
	
	#pragma omp parallel num_threads(3) private(tid)
	{
		for (int j=0; j<3; j++){
		
			int tid = omp_get_thread_num();
			cout << tid << endl;
			#pragma omp for 
			for (int i=0; i<10; i++){
				tid = omp_get_thread_num();
				cout << "i: " << i << "\t tid: " << tid << endl;
			}
		}	
	}
	
	return 0;
}
