#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "nrutil.h"
#include "math_lib.h"
#ifdef GSL
#include "GSL_tomfoolery.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#endif

void NRmat_2_GSLmat(double ** matrix, int arl, int arh, int acl, int ach, gsl_matrix * GSL_mat)
{
  for ( int i = arl; i <= arh; i++)
    {
      for ( int j = acl; j <= ach; j++)
      {
        gsl_matrix_set(GSL_mat, i-arl, j-acl, matrix[i][j]);
      }
    }
}

void GSLmat_2_NRmat(double ** matrix, int arl, int arh, int acl, int ach, gsl_matrix * GSL_mat)
{
  for ( int i = arl; i <= arh; i++)
    {
      for ( int j = acl; j <= ach; j++)
      {
        matrix[i][j] = gsl_matrix_get(GSL_mat, i-arl, j-acl);
      }
    }
}

void GSLvec_2_NRvec(double * vec, int arl, int arh, gsl_vector * GSL_vec)
{
  for ( int i = arl; i <= arh; i++)
  {
    vec[i] = gsl_vector_get(GSL_vec, i-arl);
  }
}

void NRvec_2_GSLvec(double * vec, int arl, int arh, gsl_vector * GSL_vec)
{
  for ( int i = arl; i <= arh; i++)
  {
    gsl_vector_set(GSL_vec, i-arl, vec[i]);
  }
}


void GSL_mat2vec(gsl_matrix * mat, gsl_vector * vec)
{
  int i, j;
  for ( j = 0; j < mat->size2; j++)
  {
    for ( i = 0; i < mat->size1; i++)
    {
      gsl_vector_set(vec, j*(mat->size1) + i, gsl_matrix_get(mat, i, j) );
    }
  }
}

void GSL_vec2mat(gsl_matrix * mat, gsl_vector * vec)
{
  int i, j;
  for ( j = 0; j < mat->size2; j++)
  {
    for ( i = 0; i < mat->size1; i++)
    {
      gsl_matrix_set(mat, i, j, gsl_vector_get(vec, j*(mat->size1) + i) );
    }
  }
}


void G_func_GSL(gsl_matrix * A, gsl_matrix * X, gsl_vector * b, gsl_vector * b_work, gsl_vector * x_vec)
{
  int ret_int;
  double alpha, beta;
  GSL_mat2vec(X, x_vec);

  gsl_vector_memcpy(b_work, b);
  alpha = 1;
  beta = -1;
  ret_int = gsl_blas_dgemv(CblasNoTrans, alpha, A, x_vec, beta, b_work);

  beta = 0;
  ret_int = gsl_blas_dgemv(CblasTrans, alpha, A, b_work, beta, x_vec);

}

void matrix_shrink_GSL(gsl_matrix * X, gsl_matrix * U, gsl_vector * S, gsl_matrix * V, gsl_vector * work, gsl_matrix * S_work, gsl_matrix * nn_work, double nu)
{
  int ret_int, i ;
  double alpha, beta;
  gsl_matrix_memcpy(U, X);
  ret_int = gsl_linalg_SV_decomp(U, V, S, work);

  vector_shrink_GSL(S, nu);
  for ( i = 0; i < S->size; i++)
  {
    gsl_matrix_set(S_work, i, i, gsl_vector_get(S, i) );
    // gsl_matrix_set(S_work, i, i, gsl_vector_get(S, i)/gsl_vector_get(S, 0) );

  }
  alpha = 1;
  beta = 0;
  ret_int = gsl_blas_dgemm(CblasNoTrans, CblasTrans, alpha, S_work, V, beta, nn_work);
  ret_int = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, U, nn_work, beta, X);
}

void vector_shrink_GSL(gsl_vector * S, double nu)
{
  int i ;
  for ( i = 0; i < S->size; i++)
  {
    if (gsl_vector_get(S, i) - nu > 0 )
    {
      gsl_vector_set(S, i, gsl_vector_get(S, i) - nu);
    }
    else
    {
      gsl_vector_set(S, i, 0);
    }

  }
}
