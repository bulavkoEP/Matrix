#ifndef MATR
#include "Matrix.h"
#define MATR
#endif

#include <iostream>

using namespace std;

Matrix::Matrix(int n, int m) {
    this->n = n; this->m = m;
    mat = new double[n * m];
}

Matrix::Matrix() {
    this->n = 0; this->m = 0; mat = nullptr;
}

double Matrix::get(int i, int j) {
    if (i > n || j > m) return -1;
    return mat[i * m + j];
}

int Matrix::set(int i, int j, double val) {
    if (i > n || j > m) return -1;
    mat[i * m + j] = val;
    return 1;
}

Matrix* Matrix::mult_by(Matrix* b) {
    if (b->n != m) return new Matrix(0, 0);
    Matrix* res = new Matrix(n, b->m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < b->m; ++j) {
            double c = 0;
            for (int r = 0; r < m; ++r) {
                c += get(i, r) * b->get(r, j);
            }
            res->set(i, j, c);
        }
    }
    return res;
}

Matrix* Matrix::substract(Matrix* b) {
    if (b->n != n || b->m != m) return new Matrix(0, 0);
    Matrix* res = new Matrix(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            res->set(i, j, get(i, j) - b->get(i, j));
        }
    }
    return res;
}

double Matrix::norm() {
    double max = 0, sum = 0;
    for (int j = 0; j < m; ++j) {
        sum = 0;
        for (int i = 0; i < m; ++i) {
            sum += abs(get(i, j));
        }
        if (sum > max) max = sum;
    }
    return max;
}