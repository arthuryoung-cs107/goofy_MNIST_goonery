#ifndef math_lib_H  /* Include guard */
#define math_lib_H

// auxiliary functions
void read_csv_matrix(char filename[], double ** mat_out , int m, int n );
void write_csv_matrix(char filename[], double ** mat, int m, int n);

// relevant math functions
void zerom_init(double **a, int arl, int arh, int acl, int ach);
void zeromint_init(int **a, int arl, int arh, int acl, int ach);
void zerov_init(double *a, int arl, int arh);
void zerovint_init(int *a, int arl, int arh);
void mat2vec(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh); // row stacking
void vec2mat(double * a, int arl, int arh, double ** b, int brl, int brh, int bcl, int bch); // row stacking
void d_mat_vec_mult(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c);
void dvector_add(double * a, int arl, int arh, double * b, int sign, double * result); //assumes equivalent indexing regime
void dv_scalmult(double * a, int arl, int arh, double scalar, double * result); // potentially more fitting as a scalar?
void dm_scalmult(double ** a, int arl, int arh, int acl, int ach, double scalar, double ** result);
void dmatrix_mult(double ** a, int arl, int arh, int acl, int ach, double ** b, int brl, int brh, int bcl, int bch, double ** c);
double dot(double * a, int arl, int arh, double * b, int brl, int brh);


#endif // math_lib_H
