#include "Cholesky.h"

#define EPS 1e-15

int get_inverse(Matrix* mat, Matrix* res, Matrix* decomp, double* y, double* x) {
    get_cholesky_decomp(mat, decomp);
    double norm = mat->norm();

    for (int i = 0; i < mat->n; ++i) {
        for (int j = 0; j < mat->n; ++j) {
            x[j] = 0; y[j] = 0;
        }

        for (int j = 0; j < mat->n; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                //sum += decomp->get(j, k) * y[k];
                sum += decomp->mat[j * decomp->m + k] * y[k];
            }
            double b = 0;
            if (j == i) b = 1;
            if (abs(norm) < 1e-20 || abs(decomp->mat[j * decomp->m + j]) / norm < EPS) {
                return -1;
            }
            //y[j] = 1.0 / decomp->get(j, j) * (b - sum);
            y[j] = 1.0 / decomp->mat[j * decomp->m + j] * (b - sum);
        }

        for (int j = 0; j < mat->n; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                //sum += decomp->get(mat->n - k - 1, mat->n - j - 1) * x[k];
                sum += decomp->mat[(mat->n - k - 1) * decomp->m + mat->n - j - 1] * x[k];
            }
            //x[j] = 1.0 / decomp->get(mat->n - j - 1, mat->n - j - 1) * (y[mat->n - j - 1] - sum);
            x[j] = 1.0 / decomp->mat[(mat->n - j - 1) * decomp->m + mat->n - j - 1] * (y[mat->n - j - 1] - sum);
        }
        for (int j = 0; j < mat->n; ++j) {
            //res->set(j, i, x[mat->n - j - 1]);
            res->mat[j * res->m + i] = x[mat->n - j - 1];
        }
    }
    return 1;
}

void get_cholesky_decomp(Matrix* mat, Matrix* decomp) {
    for (int i = 0; i < mat->n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++){
                //s += decomp->get(i, k) * decomp->get(j, k);
                s += decomp->mat[i * decomp->m + k] * decomp->mat[j * decomp->m + k];
            }
            
            //if (i == j) decomp->set(i, j, sqrt(mat->get(i, i) - s));
            if (i == j) decomp->mat[i * decomp->m + j] = sqrt(mat->mat[i * mat->m + i] - s);
            //else decomp->set(i, j, 1.0 / decomp->get(j, j) * (mat->get(i, j) - s));
            else decomp->mat[i * decomp->m + j] = 1.0 / decomp->mat[j * decomp->m + j] * (mat->mat[i * mat->m + j] - s);
        }
}