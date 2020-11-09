#ifndef MATR
#include "Matrix.h"
#define MATR
#endif

#include <math.h>
#include <iostream>

using namespace std;

int get_inverse(Matrix* mat, Matrix* res, Matrix* decomp, double* x, double* y);
void get_cholesky_decomp(Matrix* mat, Matrix* decomp);