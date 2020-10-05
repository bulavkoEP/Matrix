#ifndef MATR
#include "Matrix.h"
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

    Matrix* matrix = new Matrix(n, n);
    if (k != 0) {
        generate_matrix_from_formula(matrix, k);
    } else {
        read_matrix_from_file(matrix, filename);
    }
    cout << endl << matrix->mat[5 * 2 + 2] << endl;

    Matrix* res = new Matrix(n, n);
    Matrix* decomp = new Matrix(n, n);
    double* x = new double[n];
    double* y = new double[n];
    clock_t t_start = clock();
    //algo
    get_inverse(matrix, res, decomp, x, y);
    cout << "TIME: " << (double) (clock() - t_start) << endl;

    cout << "MATRIX: " << endl;
    print(matrix, m, cout);
    cout << endl;
    cout << "DECOMP: " << endl;
    print(decomp, m, cout);
    cout << endl;

    cout << "RES: " << endl;
    print(res, m, cout);
    cout << endl;

    Matrix* E = matrix->E(n);
    Matrix* mult = matrix->mult_by(res);
    Matrix* substracted = mult->substract(E);
    cout << "MULT: " << endl;
    print(mult, m, cout);
    cout << endl;
    cout << "NORM: " << substracted->norm() << endl;

    delete matrix; delete decomp; delete res; delete E;
    delete substracted; delete[] x; delete[] y; delete mult;
}