// Author: David Murphy (dmurphy@phys.columbia.edu)
// Modified: 02/19/2015

#include<cstdio>
#include<cstdlib>
#include<vector>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_eigen.h>

// Must use "gcc -lgsl -lgslcblas foo.c -o foo" to compile with gcc

class Matrix
{
  // Every instance of a matrix has this data associated with it
  const int nrows;
  const int ncols;
  double* data;

public:
  // default constructor (nrows = ncols = 0, data is empty)
  Matrix() : nrows(), ncols(), data() {};

  //Constructor for MxN matrix of zeros
  Matrix(int M, int N) : nrows(M), ncols(M), { data = new double[M*N]; };

  // Constructor for MxN matrix copied from a GSL matrix object
  Matrix(gsl_matrix* M) : nrows(M->size1), ncols(M->size2)
  {
	data = new double[nrows*ncols]
  }
}
