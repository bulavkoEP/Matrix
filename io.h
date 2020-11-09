#ifndef MATR
#include "Matrix.h"
#define MATR
#endif
#include <iostream>
#include <fstream>

double fun(int k, int n, int i, int j);

int read_matrix_from_file(Matrix* matrix, string filename);

void generate_matrix_from_formula(Matrix* matrix, int k);

void print(Matrix* matrix, int m, ostream& of);

double fun(int k, int n, int i, int j);