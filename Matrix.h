#include <string>
#include <iostream>
using namespace std;

class Matrix {
    public:
    int n, m;
    double* mat;
    Matrix(int n, int m);
    Matrix();
    double get(int i, int j);
    int set(int i, int j, double val);
    ~Matrix() {
        delete[] mat;
    }

    Matrix* E(int n) {
        Matrix* res = new Matrix(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res->set(i, j, 0);
            }
            res->set(i, i, 1);
        }
        return res;
    }
    Matrix* mult_by(Matrix* b);
    Matrix* substract(Matrix* b);
    double norm();
    private:
};