#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "aux_functions.h"
#include "nrutil.h"

// auxiliary functions for problem

void read_csv_matrix(char filename[], double ** mat_out , int m, int n )
{
  int row_count, col_count, i, j;
  // double ** mat_out = dmatrix(0, m-1, 0, n-1);
  char extract[1000];
  char *rest;
  char *rest2;
  double check;
  zerom_init(mat_out, 0, m-1, 0, n-1);

  FILE * file_stream = fopen(filename, "r");
  row_count = 0;
  while (fgets(extract, sizeof (extract), file_stream))
  {
    col_count = 0;
    if (row_count != 0)
    {
      check = strtod(extract, &rest);
      while (col_count < n)
      {
        mat_out[row_count-1][col_count] = strtod(rest + 1, &rest2);
        memcpy(rest, rest2, sizeof(extract));
        col_count++;
      }
    }
    row_count++;
  }
  fclose(file_stream);
}

void write_csv_matrix(char filename[], double ** mat, int m, int n)
{
  int i, j;

  FILE * out_file  = fopen(filename, "w");
  for ( i = 0; i < n; i++)
  {
    fprintf(out_file, ",%d", i);
  }
  fprintf(out_file, "\n");

  for ( i = 0; i < m; i++)
  {
    fprintf(out_file, "%d", i);
    for ( j = 0; j < n; j++)
    {
      fprintf(out_file, ",%f", mat[i][j]);
    }
    fprintf(out_file, "\n");
  }
  fclose(out_file);
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
