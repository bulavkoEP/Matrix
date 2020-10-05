#include "Cholesky.h"

void get_inverse(Matrix* mat, Matrix* res, Matrix* decomp, double* y, double* x) {
    get_cholesky_decomp(mat, decomp);
    
    for (int i = 0; i < mat->n; ++i) {
        for (int j = 0; j < mat->n; ++j) {
            x[j] = 0; y[j] = 0;
        }

        for (int j = 0; j < mat->n; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += decomp->get(j, k) * y[k];
            }
            double b = 0;
            if (j == i) b = 1;
            y[j] = 1.0 / decomp->get(j, j) * (b - sum);
        }
        /*
        for (int j = mat->n - 1; j >= 0; --j) {
            double sum = 0;
            for (int k = j + 1; k < mat->n; ++k) {
                sum += decomp->get(k, j) * x[k];
            }
            x[j] = 1.0 / decomp->get(j, j) * (y[j] - sum);
        }*/

        for (int j = 0; j < mat->n; ++j) {
            double sum = 0;
            for (int k = 0; k < j; ++k) {
                sum += decomp->get(mat->n - k - 1, mat->n - j - 1) * x[k];
            }
            x[j] = 1.0 / decomp->get(mat->n - j - 1, mat->n - j - 1) * (y[mat->n - j - 1] - sum);
        }
        for (int j = 0; j < mat->n; ++j) {
            res->set(j, i, x[mat->n - j - 1]);
        }
    }
}

void get_cholesky_decomp(Matrix* mat, Matrix* decomp) {
    for (int i = 0; i < mat->n; i++)
        for (int j = 0; j < (i+1); j++) {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += decomp->get(i, k) * decomp->get(j, k);
            
            if (i == j) decomp->set(i, j, sqrt(mat->get(i, i) - s));
            else decomp->set(i, j, 1.0 / decomp->get(j, j) * (mat->get(i, j) - s));
        }
}