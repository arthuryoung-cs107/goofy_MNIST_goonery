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

int main()
{
  /////////////////////////////////////////////
  ///////////// PROBLEM 2 solver //////////////
  /////////////////////////////////////////////

  char pierce_filename[50];
  char specfile[300];
  snprintf(pierce_filename, 49, "../pierce.txt");

  int i, j;
  int m = 100;
  int n = 200;
  int C_vec_ind, G_vec_ind, M_vec_ind;
  int ** pierce_mat = Pierce_mat_gen(pierce_filename, m, n);
  double ** G = G_mat_gen(pierce_mat, 1, m, 1, n);

  // problem parameters
  double p0 = 10;
  double h = 36.6e-2;
  double c = 3.43e2;
  double w = 100*M_PI;
  double del_T_mult = 1;
  // double del_T = del_T_mult*h/(2*c);
  double del_T = 5e-5;
  double t_it = del_T;
  double t_max = 1.01;
  double scale_update = del_T*del_T*c*c/(h*h);
  double hear_tol = 1e-3;
  int n_step = 1;
  int * C_ind = ivector(1, 2);
  int * G_ind = ivector(1, 2);
  int * M_ind = ivector(1, 2);
  double * P_0_vec = dvector(1, m*n);
  double * P_1_vec = dvector(1, m*n);
  double * P_k_vec = dvector(1, m*n);
  double * P_k_old = dvector(1, m*n);
  double * P_k_new = dvector(1, m*n);
  double * time_vec = dvector(1, 4);
  double * Pvec_scratch = dvector(1, m*n );
  double ** Pmat_scratch = dmatrix(1, m, 1, n);
  double ** Pmat_out = dmatrix(1, m, 1, n);
  int ** S_ind = imatrix(1, 2, 1, 2);
  int time_count = 1;
  int C_switch = 0;
  int G_switch = 0;
  int M_switch = 0;

  C_ind[1] = 35 + 1;
  C_ind[2] = 73 + 1;
  G_ind[1] = 61 + 1;
  G_ind[2] = 109 + 1;
  M_ind[1] = 91 + 1;
  M_ind[2] = 188 + 1;
  C_vec_ind = n*(C_ind[1] - 1) + C_ind[2];
  G_vec_ind = n*(G_ind[1] - 1) + G_ind[2];
  M_vec_ind = n*(M_ind[1] - 1) + M_ind[2];
  time_vec[1] = 1.5e-2;
  time_vec[2] = 1.05e-1;
  time_vec[3] = 5.05e-1;
  time_vec[4] = 1.005;
  S_ind[1][1] = 58;
  S_ind[1][2] = 61;
  S_ind[2][1] = 16;
  S_ind[2][2] = 19;
  zerov_init(P_0_vec, 1, m*n);
  zerov_init(P_1_vec, 1, m*n);
  S_solve(S_ind, P_1_vec, Pmat_scratch, m, n, p0, w, t_it);

  dvector_cpy(P_1_vec, 1, m*n, P_k_vec, 1, m*n);
  dvector_cpy(P_0_vec, 1, m*n, P_k_old, 1, m*n);

  printf("delta t: %f\n", del_T);

  memset(specfile, 0, 299);
  snprintf(specfile, 300, "Hear_times.dat");
  FILE * hear_time_file = fopen(specfile, "w");

  memset(specfile, 0, 299);
  snprintf(specfile, 300, "position_pressures_over_time.dat");
  FILE * position_pressure_file = fopen(specfile, "w");

  while (t_it < t_max)
  {
    n_step++;
    t_it = t_it + del_T;
    // d_mat_vec_mult(G, 1, m*n, 1, m*n, P_k_vec, 1, m*n, Pvec_scratch); // perform spacial update
    // d_mat_vec_mult_special(G, 1, m*n, 1, m*n, P_k_vec, 1, m*n, Pvec_scratch, m, n); // special matrix multiplication
    d_mat_vec_mult_special2(G, 1, m*n, 1, m*n, P_k_vec, 1, m*n, Pvec_scratch, m, n); // special matrix multiplication
    dv_scalmult(Pvec_scratch, 1, m*n, scale_update, Pvec_scratch);
    dvector_add(Pvec_scratch, 1, m*n, P_k_vec, 2, Pvec_scratch);
    dvector_add(Pvec_scratch, 1, m*n, P_k_old, -1, Pvec_scratch);
    clean_up(S_ind, pierce_mat, Pvec_scratch, m, n, Pmat_scratch, P_k_new, p0, w, t_it);

    dvector_cpy(P_k_vec, 1, m*n, P_k_old, 1, m*n);
    dvector_cpy(P_k_new, 1, m*n, P_k_vec, 1, m*n);
    if (n_step%50 == 0)
    {
      printf("t = %f, t_crit = %f \n", t_it, time_vec[time_count]);
      fprintf(position_pressure_file, "%e %e %e %e \n ", t_it, P_k_vec[C_vec_ind], P_k_vec[G_vec_ind], P_k_vec[M_vec_ind]);
    }
    if (1000*fabs(t_it - time_vec[time_count]) < del_T)
    {
      vec2mat(P_k_vec, 1, m*n, Pmat_out, 1, m, 1, n);
      memset(specfile, 0, 299);
      snprintf(specfile, 300, "Pressure_field_time%d.dat", time_count);
      FILE * Pmat_out_file = fopen(specfile, "w");
      for ( i = 1; i <= m; i++)
      {
        for ( j = 1; j <= n; j++)
        {
          fprintf(Pmat_out_file, "%e ", Pmat_out[i][j]);
        }
        fprintf(Pmat_out_file, "\n");
      }
      fclose(Pmat_out_file);
      time_count++;
    }
    if ( (fabs(P_k_vec[C_vec_ind]) >= hear_tol) && (C_switch == 0))
    {
      fprintf(hear_time_file, "C: %e \n", t_it);
      C_switch = 1;
    }
    if ( (fabs(P_k_vec[G_vec_ind]) >= hear_tol) && (G_switch == 0 ) )
    {
      fprintf(hear_time_file, "G: %e \n", t_it);
      G_switch = 1;
    }
    if  ( (fabs(P_k_vec[M_vec_ind]) >= hear_tol) && (M_switch == 0) )
    {
      fprintf(hear_time_file, "M: %e \n", t_it);
      M_switch = 1;
    }

  }
  fclose(hear_time_file);
  fclose(position_pressure_file);


  return 0;
}
