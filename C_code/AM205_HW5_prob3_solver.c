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
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#endif

int main()
{
  #ifdef GSL
    printf("Employing GSL matrix algebra packages\n");
  #endif
  int i, j, m;
  int count = 0;
  int out_check;
  int signum;
  int n = 101;
  int n_2 = n- 2;
  double norm_it1, norm_it2;
  double tol = 1e-10;
  double err = 1;
  double min = -1;
  double max = 1;
  double del = (max-min)/( (double)(n - 1) );
  double x_it = 0;
  double h_inv_sqr = pow(del, -2);

  double * x_vec = dvector(1, n);

  // printf("del: %f\n", del);
  for ( i = 1 ; i <= n; i++)
  {
    x_vec[i] = x_it + min;
    x_it = x_it + del;
  }

  gsl_matrix * Jacobian_it = gsl_matrix_alloc(n_2, n_2);
  gsl_matrix * LU_it = gsl_matrix_alloc(n_2, n_2);
  gsl_permutation * perm_it = gsl_permutation_alloc(n_2);
  gsl_vector * u_it = gsl_vector_alloc(n_2);
  gsl_vector * f_it = gsl_vector_alloc(n_2);
  gsl_vector * u_old = gsl_vector_alloc(n_2);
  gsl_vector * x_solve = gsl_vector_alloc(n_2);
  gsl_vector * x_diff = gsl_vector_alloc(n_2);
  gsl_vector_set_zero(u_it);
  gsl_vector_set_zero(u_old);

  while(err > tol)
  {
    Jacobian_gen_GSL(Jacobian_it, u_it, n_2, h_inv_sqr);
    F_gen_GSL(f_it, u_it, n_2, h_inv_sqr);
    gsl_blas_dscal(-1, f_it);
    gsl_matrix_memcpy(LU_it, Jacobian_it);
    gsl_linalg_LU_decomp(LU_it, perm_it, &signum);
    gsl_linalg_LU_solve(LU_it, perm_it, f_it, x_solve);
    gsl_vector_add(u_it, x_solve);
    gsl_vector_memcpy(x_diff, u_it);
    gsl_vector_sub(x_diff, u_old);

    err = (gsl_blas_dnrm2(x_diff) )/(gsl_blas_dnrm2(u_old) );

    gsl_vector_memcpy(u_old, u_it);
    count++;
  }

  printf("error: %e\n", err);
  printf("count: %d\n", count);
  char specfile[300];

  memset(specfile, 0, 299);
  snprintf(specfile, 300, "u_converged.dat");
  FILE * u_out_file  = fopen(specfile, "w");
  for ( i = 0; i < n_2; i++)
  {
    fprintf(u_out_file, "%e %e\n", x_vec[i + 2], gsl_vector_get(u_it, i) );
  }
  fclose(u_out_file);



  return 0;
}
