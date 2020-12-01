#include "Cholesky.h"

#define EPS 1e-15


int get_inverse(double* mat, double* res, double* d, int n) {
    if (cholesky_decomp(mat, res, d, n) == -1) return -1;
    
    double s = 0;
    
    for (int i = 0; i < n; ++i)
		for (int j = i; j >= 0; --j)
		{
			s = (double)(i == j);
			for (int k = j + 1; k <= i; ++k)
				s -= mat[ind(i, k)] * res[ind(j, k)];
			mat[ind(i, j)] = s / res[ind(j, j)];
		}
	
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
		{
			s = 0.0;
			for (int k = (i > j) ? i : j; k < n; ++k)
				s += d[k] * mat[ind(k, j)] * mat[ind(k, i)];
			res[ind(i, j)] = s;
		}
	
    return 1;
}

int cholesky_decomp(double* mat, double* res, double* d, int n) {
    double s = 0;
    double no = norm(mat, n);
    for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
		{
			s = mat[ind(j, i)];
			for (int k = 0; k < i; ++k)
				s -= res[ind(i, k)] * res[ind(j, k)] * d[k];

			if (i == j)
			{
				if (s < 0)
				{
					s = -s;
					d[i] = -1.0;
				}
				else
					d[i] = 1.0;

				if (fabs(s / no) < EPS)
					return -1;

				s = sqrt(s);
				res[ind(i, i)] = s;
			}
			else
				res[ind(j, i)] = s / (res[ind(i, i)] * d[i]);
		}
    return 1;
}