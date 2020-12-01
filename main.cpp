#ifndef MATR
#include "functions.h"
#define MATR
#endif
#include "io.h"
#include <time.h>

#include "Cholesky.h"

using namespace std;

int main(int argc, char * argv[]) {
    string filename = "";
    int n, m, k;
    if (argc == 5) {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
        filename = argv[4];
    } else if (argc == 4) {
        n = stoi(argv[1]);
        m = stoi(argv[2]);
        k = stoi(argv[3]);
    } else {
        cout << "Usage: ./main n m k filename or ./main n m k when k != 0";
        return 0;
    }

    double* matrix = new double[n * (n + 1) / 2];
    if (k != 0) {
        generate_matrix_from_formula(matrix, n, k);
    } else {
        if (read_matrix_from_file(matrix, n, filename) != 1) {
            delete[] matrix;
            return 0;
        }
    }

    cout << "MATRIX: " << endl;
    print(matrix, n, m, cout);
    cout << endl;

    double* res = new double[n * (n + 1) / 2];
    double* d = new double[n];
    clock_t t_start = clock();

    int r = get_inverse(matrix, res, d, n);

    if (r == -1) {
        cout << "det is zero" << endl;
         delete[] matrix; delete[] res;
         delete[] d; 
         return 0;
    }

    cout << "TIME: " << (double) (clock() - t_start) / CLOCKS_PER_SEC << endl;

    cout << "A: " << endl;
    print(matrix, n, m, cout);
    cout << "RES: " << endl;
    print(res, n, m, cout);
    cout << endl;

    if (k != 0) {
        generate_matrix_from_formula(matrix, n, k);
    } else {
        if (read_matrix_from_file(matrix, n, filename) != 1) {
            delete matrix;
            return 0;
        }
    }

    cout << "NORM: " << mult_err(matrix, res, n) << endl;

    delete[] matrix; delete[] res; delete[] d; 
}