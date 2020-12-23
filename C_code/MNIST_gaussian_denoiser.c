#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "nrutil.h"
#include "aux_functions.h"
#ifdef GSL
#include "GSL_tomfoolery.h"
#include "line_search.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

int main()
{
  #ifdef GSL
    printf("Employing GSL matrix algebra packages\n");
  #endif

  char mnist_filename[50];
  char specfile[300];

  int i, j, ret_int, rank, p, it, mu_count;
  int m = 28;
  int n = 28;
  int max_it = 500;
  int mnist_max = 1;

  double A_std, mu_high, mu, nu, cond, tau_max, p_d, A_S1;
  double m_d = m;
  double n_d = n;
  double S00_1;
  double alpha;
  double beta;
  double noise_var = 0.01;
  double rank_tol = 1e-10;

  double ** mnist_double_matrix = dmatrix(0, m-1, 0, n-1);

  gsl_matrix * X00 = gsl_matrix_alloc(m, n);
  gsl_matrix * U00 = gsl_matrix_alloc(m, n);
  gsl_vector * S00 = gsl_vector_alloc(n);
  gsl_matrix * V00 = gsl_matrix_alloc(n, n);
  gsl_matrix * V00_T = gsl_matrix_alloc(n, n);
  gsl_matrix * X = gsl_matrix_alloc(m, n);
  gsl_matrix * U = gsl_matrix_alloc(m, n);
  gsl_vector * S = gsl_vector_alloc(n);
  gsl_matrix * V = gsl_matrix_alloc(n, n);
  gsl_matrix * V_T = gsl_matrix_alloc(n, n);

  gsl_matrix * Xk = gsl_matrix_alloc(m, n);
  gsl_matrix * Uk = gsl_matrix_alloc(m, n);
  gsl_vector * Sk = gsl_vector_alloc(n);
  gsl_matrix * Vk = gsl_matrix_alloc(n, n);
  gsl_matrix * Vk_T = gsl_matrix_alloc(n, n);
  gsl_matrix * Xk_work = gsl_matrix_alloc(m, n);
  gsl_matrix * Xk_work2 = gsl_matrix_alloc(m, n);
  gsl_matrix * Xk_last = gsl_matrix_alloc(m, n);
  gsl_matrix * S_mat_work = gsl_matrix_alloc(n, n);
  gsl_matrix * nn_work = gsl_matrix_alloc(n, n);


  gsl_vector * work = gsl_vector_alloc(n);
  gsl_vector * Xvec = gsl_vector_alloc(m*n);
  gsl_vector * Xvec_work = gsl_vector_alloc(m*n);
  gsl_vector * Xvec_work2 = gsl_vector_alloc(m*n);

  gsl_rng * T1 = gsl_rng_alloc(gsl_rng_taus);

  for ( int file_it = 0; file_it < mnist_max; file_it++)
  {
    memset(mnist_filename, 0, 49);
    snprintf(mnist_filename, 49, "../MNist_csv_raw/MNIST%d.csv", file_it);
    read_csv_matrix(mnist_filename, mnist_double_matrix, m, n);
    NRmat_2_GSLmat(mnist_double_matrix, 0, m-1, 0, n-1, X00);

    gsl_matrix_memcpy(U00, X00);

    ret_int = gsl_linalg_SV_decomp(U00, V00, S00, work);
    ret_int = gsl_matrix_transpose_memcpy(V00_T, V00);
    S00_1 = gsl_vector_get(S00, 0);

    rank = 0;
    for ( i = 0; i < n; i++)
    {
      if (gsl_vector_get(S00, i) > rank_tol)
      {
        rank++;
      }

    }
    printf("MNIST%d native rank: %d out of %d \n", file_it,rank, n);
    p = 3*rank*(m + n - rank);
    p_d = p;
    A_std = sqrt(1/(p_d));

    gsl_matrix * A = gsl_matrix_alloc(p, m*n);
    gsl_matrix * A_U = gsl_matrix_alloc(p, m*n);
    gsl_vector * A_s = gsl_vector_alloc(m*n);
    gsl_matrix * A_V = gsl_matrix_alloc(m*n, m*n);
    gsl_vector * A_work = gsl_vector_alloc(m*n);

    gsl_vector * b = gsl_vector_alloc(p);
    gsl_vector * b_k = gsl_vector_alloc(p);
    gsl_vector * b_work = gsl_vector_alloc(p);

    for ( i = 0; i < p; i++)
    {
      for ( j = 0; j < m*n; j++)
      {
        gsl_matrix_set(A, i, j, gsl_ran_gaussian(T1, A_std));
      }
    }

    tau_max = 2/((1 + sqrt(m_d*n_d/p_d) )* (1 + sqrt(m_d*n_d/p_d) ) );
    printf("tau_max %f\n", tau_max);

    for ( j = 0; j < n; j++) // down columns
    {
      for ( i = 0; i < m; i++) // down rows
      {
        gsl_vector_set(Xvec, j*(m) + i, gsl_matrix_get(X00, i, j)/S00_1 );
        gsl_matrix_set(X, i, j, gsl_matrix_get(X00, i, j)/S00_1);
      }
    }

    alpha = 1;
    beta = 0;
    ret_int = gsl_blas_dgemv(CblasNoTrans, alpha, A, Xvec, beta, b);
    ret_int = gsl_blas_dgemv(CblasTrans, alpha, A, b, beta, Xvec_work);
    gsl_vector_memcpy(b_k, b);
    gsl_matrix_set_zero(S_mat_work);

    int breg_max = 4;
    int max_it = 1000;
    double eta = 0.75;
    double mu_low = 1e-4;
    double gtol = 1e-8;
    double xtol = 1e-10;

    mu = eta*gsl_blas_dnrm2(Xvec_work);
    gsl_matrix_set_zero(Xk);

    FILE * interim_spec  = fopen("./MNIST_interim/interim_spec_file.dat", "w");
    for (int breg_it = 1; breg_it < breg_max; breg_it++)
    {
      mu_count = 0;
      while (mu > mu_low) // fixed point continuation
      {
        mu_count++;
        cond = 1;
        it = 0;
        nu = mu*tau;
        while ((cond > xtol)&&(it < max_it)) // fixed point iteration
        {
          it++;
          gsl_matrix_memcpy(Xk_last, Xk);

          G_func_GSL(A, Xk, b_k, b_work, Xvec);
          GSL_vec2mat(Xk_work, Xvec);

          tau = line_search_tau(A, Xk, b_k, Xvec, Xvec_work, Xvec_work2, b_work, 2*tau_max);
          nu = tau*mu;


          gsl_matrix_scale(Xk_work, tau );
          gsl_matrix_sub(Xk, Xk_work);

          matrix_shrink_GSL(Xk, Uk, Sk, Vk, work, S_mat_work, nn_work, nu);
          gsl_matrix_memcpy(Xk_work, Xk);
          gsl_matrix_sub(Xk_work, Xk_last);
          GSL_mat2vec(Xk_work, Xvec);
          GSL_mat2vec(Xk_last, Xvec_work);
          cond = gsl_blas_dnrm2(Xvec_work);
          if (1 > cond)
          {
            cond = 1;
          }
          cond = gsl_blas_dnrm2(Xvec) / cond;
        }
        printf("breg_it = %d , mu_it = %d, mu = %f, it_final = %d, cond = %f \n", breg_it, mu_count, mu, it, cond);

        mu = mu*eta;
        gsl_matrix_memcpy(Xk_work, Xk);
        gsl_matrix_scale(Xk_work, (double) S00_1 );
        GSLmat_2_NRmat(mnist_double_matrix, 0, m-1, 0, n-1, Xk_work);

        memset(mnist_filename, 0, 49);
        snprintf(mnist_filename, 49, "./MNIST_interim/MNIST%d_CRR_%d_%d.csv", file_it, breg_it, mu_count);

        write_csv_matrix(mnist_filename, mnist_double_matrix, m, n);

      }
      fprintf(interim_spec, "%d\n", mu_count);

      GSL_mat2vec(Xk, Xvec);
      alpha = 1;
      beta = 0;
      ret_int = gsl_blas_dgemv(CblasNoTrans, alpha, A, Xvec, beta, b_work);
      ret_int = gsl_vector_sub(b_k, b_work);
      ret_int = gsl_vector_add(b_k, b);
      ret_int = gsl_blas_dgemv(CblasTrans, alpha, A, b_k, beta, Xvec_work);
      mu = eta*gsl_blas_dnrm2(Xvec_work);
    }
    fclose(interim_spec);
    gsl_matrix_scale(Xk, (double) S00_1 );
    GSLmat_2_NRmat(mnist_double_matrix, 0, m-1, 0, n-1, Xk);

    memset(mnist_filename, 0, 49);
    snprintf(mnist_filename, 49, "../MNist_csv_CRR/MNIST%d_CRR.csv", file_it);

    write_csv_matrix(mnist_filename, mnist_double_matrix, m, n);

    gsl_vector_free(b);
    gsl_vector_free(b_k);
    gsl_vector_free(b_work);
    gsl_matrix_free(A);
  }



  return 0;
}
