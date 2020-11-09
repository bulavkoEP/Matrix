#include "io.h"


int read_matrix_from_file(Matrix* matrix, string filename) {
    ifstream in;
    in.open(filename);
    if (!in) {
        cout << "Error opening file" << endl;
        return -1;
    }
    double tmp;
    
    for (int i = 0; i < matrix->n; ++i) {
        for (int j = 0; j < matrix->m; ++j) {
            in >> tmp;
            matrix->set(i, j, tmp);
        }
    }
    in.close();
    return 1;
}

void generate_matrix_from_formula(Matrix* matrix, int k) {
    for (int i = 0; i < matrix->n; ++i) {
        for (int j = 0; j < matrix->n; ++j) {
            matrix->set(i, j, fun(k, matrix->n, i, j));
        }
    }
}

void print(Matrix* matrix, int m, ostream& of) {
    for (int i = 0; i < min(matrix->n, m); ++i) {
        for (int j = 0; j < min(matrix->m, m); ++j) {
            of << matrix->get(i, j) << " ";
        }
        of << endl;
    }
}

double fun(int k, int n, int i, int j) {
    switch(k) {
        case 1: {
            return n - max(i, j) + 1;
        }
        case 2: {
            return max(i, j);
        }
        case 3: {
            return abs(i - j);
        }
        case 4: {
            return 1 / (i + j - 1);
        }
        default: {
            return 0;
        }
    }
} 