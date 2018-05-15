#include <iostream>
#include <omp.h>
using namespace std;

int main(){
	int a;
	#pragma omp parallel
	//a=0;
	
	for (int i=0; i<4; i++){
		#pragma omp parallel for
		for (int j=0; j<4; j++){
			cout << j << endl;
		}
	}
	
	
	
	return 0;
}
