#include <iostream>

using namespace std;

int main(){
    int *a = new int(5);
    int size = sizeof(a)/sizeof(a[0]);
    cout << size << endl;
    return 0;
}
