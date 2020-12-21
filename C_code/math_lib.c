#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "math_lib.h"
#include "nrutil.h"

// auxiliary functions for problem

int ** Pierce_mat_gen(char pierce_filename[] , int m, int n )
{
  // note: we ignore edge and corner values, since they're outside the domain of the actual problem. Asusme
  int row_count, col_count, i, j;
  int ** pierce_mat = imatrix(1, m, 1, n);
  char extract[1000];
  zeromint_init(pierce_mat, 1, m, 1, n);

  FILE * pierce_file = fopen(pierce_filename, "r");

  while (fgets(extract, sizeof (extract), pierce_file))
  {
    row_count++;
    col_count = 0;
    for (i = 0; i < 1000; i++)
    {
      if (extract[i] == '0')
      {
        col_count++;
        pierce_mat[row_count][col_count] = 0;
        // printf("0 ");
      }
      else
      {
          if (extract[i] == '1')
          {
            col_count++;
            pierce_mat[row_count][col_count] = 1;
            // printf("1 ");
          }
      }
    }
    // printf("\n");
  }
  fclose(pierce_file);
  return pierce_mat;
}

double ** G_mat_gen(int **pierce_mat, int arl, int arh, int acl, int ach)
{
  int i, j, left_wall, right_wall, up_wall, down_wall, row_ind;
  int m = (arh-arl) +  1;
  int n = (ach-acl) +  1;
  double ** G = dmatrix(1, m*n, 1, m*n);
  zerom_init(G, 1, m*n, 1, m*n);

  // default update matrix formed, now we need to consider walls:
  for ( i = 2; i <= m-1; i++)
  {
    for ( j = 2; j <= n-1; j++)
    {
      if (pierce_mat[i][j] == 0) // only operating on non-wall spaces
      {
        row_ind =  n*(i - 1) + j;

        up_wall = 0;
        down_wall = 0;
        left_wall = 0;
        right_wall = 0;

        if (pierce_mat[i+1][j] == 1) // wall north of node
        {
          up_wall = 1;
        }
        if (pierce_mat[i-1][j] == 1) // wall south of node
        {
          down_wall = 1;
        }
        if (pierce_mat[i][j-1] == 1) // wall west of node
        {
          left_wall = 1;
        }
        if (pierce_mat[i][j+1] == 1) // wall east of node
        {
          right_wall = 1;
        }

        G[row_ind][n*(i - 1) + j] = -4 + up_wall + down_wall + left_wall + right_wall;
        G[row_ind][n*(i - 1) + j - 1] = 1 - left_wall;
        G[row_ind][n*(i - 1) + j + 1] = 1 - right_wall;
        G[row_ind][n*(i - 1) + j - n] = 1 - up_wall;
        G[row_ind][n*(i - 1) + j + n] = 1 - down_wall;
      }
    }
  }
  printf("m: %d, n: %d \n", m, n);
  return G;
}

void S_solve(int ** S_ind, double * P_vec, double ** P_mat, int m, int n, double p0, double w, double t) //overwrites vector contents with sin function
{
  int i, j;
  vec2mat(P_vec, 1, m*n, P_mat, 1, m, 1, n);
  for ( i = S_ind[1][1]; i <= S_ind[1][2]; i++)
  {
    for ( j = S_ind[2][1]; j <= S_ind[2][2]; j++)
    {
      P_mat[i][j] = p0*sin(w*t);
    }
  }
  mat2vec(P_mat, 1, m, 1, n, P_vec, 1, m*n);
}

void clean_up(int ** S_ind, int ** pierce_mat, double * P_vec, int m, int n, double ** P_mat, double * P_out, double p0, double w, double t_it)
{
  int i, j;
  // vec2mat(P_vec, 1, m*n, P_mat, 1, m, 1, n);
  // for ( i = 1; i <= m; i++)
  // {
  //   for ( j = 1; j <= n; j++)
  //   {
  //     if (pierce_mat[i][j] == 1)
  //     {
  //       P_mat[i][j] = 0;
  //     }
  //   }
  // }

  // for ( i = 1; i <= m; i++) // down rows, asserting edges to be zero for no blowup
  // {
  //   P_mat[i][1] = 0;
  //   P_mat[i][n] = 0;
  // }
  // for ( i = 1; i <= m; i++) // down cols, asserting edges to be zero for no blowup
  // {
  //   P_mat[1][i] = 0;
  //   P_mat[m][i] = 0;
  // }
  // mat2vec(P_mat, 1, m, 1, n, P_vec, 1, m*n);
  S_solve(S_ind, P_vec, P_mat, m, n, p0, w, t_it);
  dvector_cpy(P_vec, 1, m*n, P_out, 1, m*n);
}

void d_mat_vec_mult_special(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c, int m, int n) // this is a special version of matrix multiplication that I made for banded matrices such as the one in this problem. Considerably faster than the alternative.
{
  int i, j;
  double sum;

  int an = ach - acl + 1;
  int bm = brh - brl + 1;

  // int m = arh-arl + 1;
  // int n = ach-acl + 1;
  if (an != bm)
  {
    printf(" Error: matrix-vector pair of inequivalent inner-dimensions. Matrix multiplication failed \n");
  }
  else // we're handfeeding subsets of the matrix for multiplication
  {
    i = arl;
    c[i] = dot(a[i], i, i+n, b, i, i+n); // first term
    for ( i = arl+1; i <= n; i++) // first row
    {
      c[i] = dot(a[i], i-1, i+n, b, i-1, i+n);
    }
    for ( i = n + 1; i <= (m-1)*n; i++) // bulk of matrix
    {
      c[i] = dot(a[i], i-n, i+n, b, i-n, i+n);
    }
    for ( i = (n*(m-1) + 1); i <= (n*m - 1); i++) // last row
    {
      c[i] = dot(a[i], i-n, i+1, b, i-n, i+1);
    }
    i = arh;
    c[i] = dot(a[i], i-n, i, b, i-n, i);  // last element
  }
}
void d_mat_vec_mult_special2(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c, int m, int n) // this is a special version of matrix multiplication that I made for banded matrices such as the one in this problem. Considerably faster than the alternative.
{
  int i, j;
  double sum;

  int an = ach - acl + 1;
  int bm = brh - brl + 1;

  // int m = arh-arl + 1;
  // int n = ach-acl + 1;
  if (an != bm)
  {
    printf(" Error: matrix-vector pair of inequivalent inner-dimensions. Matrix multiplication failed \n");
  }
  else // we're handfeeding subsets of the matrix for multiplication
  {
    i = arl;
    // c[i] = dot(a[i], i, i+n, b, i, i+n); // first term
    c[i] = a[i][i]*b[i] + a[i][i + 1]*b[i + 1] + a[i][i + n]*b[i + n];
    for ( i = arl+1; i <= n; i++) // first row
    {
      // c[i] = dot(a[i], i-1, i+n, b, i-1, i+n);
      c[i] = a[i][i-1]*b[i-1] + a[i][i]*b[i] + a[i][i+1]*b[i+1]  + a[i][i + n]*b[i + n];
    }
    for ( i = n + 1; i <= (m-1)*n; i++) // bulk of matrix
    {
      // c[i] = dot(a[i], i-n, i+n, b, i-n, i+n);
      c[i] = a[i][i-n]*b[i-n] + a[i][i-1]*b[i-1] + a[i][i]*b[i] + a[i][i+1]*b[i+1]  + a[i][i + n]*b[i + n];
    }
    for ( i = (n*(m-1) + 1); i <= (n*m - 1); i++) // last row
    {
      // c[i] = dot(a[i], i-n, i+1, b, i-n, i+1);
      c[i] = a[i][i-n]*b[i-n] + a[i][i-1]*b[i-1] + a[i][i]*b[i] + a[i][i+1]*b[i+1];
    }
    i = arh;
    c[i] = dot(a[i], i-n, i, b, i-n, i);  // last element
    c[i] = a[i][i-n]*b[i-n] + a[i][i-1]*b[i-1] + a[i][i]*b[i] ;
  }
}




// Basic linear algebra stuff

void zerom_init(double **a, int arl, int arh, int acl, int ach)
{
  int i, j;
  for ( i = arl; i <= arh; i++)
  {
    for ( j = acl; j <= ach ; j++)
    {
      a[i][j] = 0;
    }
  }
}

void zeromint_init(int **a, int arl, int arh, int acl, int ach)
{
  int i, j;
  for ( i = arl; i <= arh; i++)
  {
    for ( j = acl; j <= ach ; j++)
    {
      a[i][j] = 0;
    }
  }
}
void zerov_init(double *a, int arl, int arh)
{
  int i;
  for ( i = arl; i <= arh; i++)
  {
    a[i] = 0;
  }
}

void zerovint_init(int *a, int arl, int arh)
{
  int i;
  for ( i = arl; i <= arh; i++)
  {
    a[i] = 0;
  }
}

void mat2vec(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh) // row stacking
{
  int rows = arh-arl + 1;
  int cols = ach-acl + 1;

  memcpy(b + brl, a[arl] + acl, rows*cols*sizeof(double));
}

void vec2mat(double * a, int arl, int arh, double ** b, int brl, int brh, int bcl, int bch) // row stacking
{
  int rows = brh-brl + 1;
  int cols = bch-bcl + 1;

  memcpy(b[brl] + bcl, a + arl, rows*cols*sizeof(double));
}

void d_mat_vec_mult(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c)
{
  int i, j;
  double sum;

  int an = ach - acl + 1;
  int bm = brh - brl + 1;

  if (an != bm)
  {
    printf(" Error: matrix-vector pair of inequivalent inner-dimensions. Matrix multiplication failed \n");
  }
  else
  {
    for ( i = arl; i <= arh; i++ ) // for each row of output vector
    {
      sum = 0;
      for(j = acl; j <= ach; j++) // each column component of i'th row of a
      {
        sum = sum + a[i][j]*b[j];
      }
      c[i] = sum;
    }
  }
}
void dvector_add(double * a, int arl, int arh, double * b, int sign, double * result) //assumes equivalent indexing regime
{
  int i;
  double sign_d = (double)sign;
  for ( i = arl; i <= arh; i++)
  {
    result[i] = a[i] +  sign_d*b[i];
  }
}

void dv_scalmult(double * a, int arl, int arh, double scalar, double * result) // potentially more fitting as a scalar?
{
  for (int i = arl; i <= arh; i++)
  {
    result[i] = a[i]*scalar;
  }
}

void dm_scalmult(double ** a, int arl, int arh, int acl, int ach, double scalar, double ** result)
{
  for (int i = arl; i <= arh; i++)
  {
    for (int j = acl; j <= ach; j++)
    {
      result[i][j] = scalar*a[i][j];
    }
  }
}

void dmatrix_mult(double ** a, int arl, int arh, int acl, int ach, double ** b, int brl, int brh, int bcl, int bch, double ** c)
{
  int i, j, k;
  double sum;
  int am = arh - arl + 1; //producing dimensions of matrices
  int an = ach - acl + 1;
  int bm = brh - brl + 1;
  int bn = bch - bcl + 1;

  if (an != bm)
  {
      printf(" Error: inner matrix dimensions inequivalent. Matrix multiplication failed. \n");
  }

  for ( i = arl; i <= arh; i++)
  {
    // printf("%s \n", "passed internal multipliation");

    for ( j = bcl; j <= bch; j++) // outer loop for resolving appropriate dimensions of output
    {
      sum = 0;
      for(k = acl; k <= ach; k++)
      {
        sum = sum + a[i][k]*b[k][j];
      }
      c[i][j] = sum;
    }
  }
}

double dot(double * a, int arl, int arh, double * b, int brl, int brh)
{
  double sum = 0;
  if ((arh - arl) != (brh - brl) )
  {
    printf("dot_product: vector dimensions inequivalent, dot product failed\n");
  }
  else
  {
    for (int i = arl; i <= arh; i++)
    {
      sum = sum + a[i]*b[i];
    }
  }
  return sum;

}
