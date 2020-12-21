#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "line_search.h"

#ifdef GSL
#include "GSL_tomfoolery.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_min.h>
#endif

struct min_func_params {gsl_matrix * A; gsl_vector * xk; gsl_vector * gk ; gsl_vector * b; gsl_vector * xvec_work; gsl_vector * b_work;};

double min_func(double tau, void * p)
{
  struct min_func_params * params = (struct min_func_params * )p;
  double norm, alpha, beta;
  int ret_int;

  gsl_vector_memcpy(params->xvec_work, params->gk );
  gsl_vector_scale( params->xvec_work , ((double) -1) * tau );
  gsl_vector_add( params->xvec_work, params->xk );
  gsl_vector_memcpy(params->b_work, params->b);

  alpha = 1;
  beta = -1;

  ret_int = gsl_blas_dgemv(CblasNoTrans, alpha, params->A, params->xvec_work, beta, params->b_work);
  norm = gsl_blas_dnrm2(params->b_work);
  return norm;
}

double line_search_tau(gsl_matrix * A, gsl_matrix * Xk, gsl_vector * bk, gsl_vector * gk, gsl_vector * xvec, gsl_vector * xvec_work, gsl_vector * b_work, double tau_max)
{
  int status;
  int iter = 0;
  int max_iter = 100;
  const gsl_min_fminimizer_type *T;
  double m = tau_max/2;
  double a = 0;
  double b = tau_max;
  double tau_return;
  gsl_min_fminimizer *s;
  GSL_mat2vec(Xk, xvec);
  gsl_function F;
  struct min_func_params F_params= {A, xvec, gk, bk, xvec_work, b_work};
  F.function = &min_func;
  F.params = &F_params;

  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, m, a, b);

  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    m = gsl_min_fminimizer_x_minimum (s);
    a = gsl_min_fminimizer_x_lower (s);
    b = gsl_min_fminimizer_x_upper (s);
    // printf("iteration = %d, m = %f \n", iter, m);

    status = gsl_min_test_interval (a, b, 0.001, 0.0);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  tau_return = m;

  gsl_min_fminimizer_free (s);

  return tau_return;
}
