#ifndef math_lib_H  /* Include guard */
#define math_lib_H

// auxiliary functions
int ** Pierce_mat_gen(char pierce_filename[] , int m, int n );
double ** G_mat_gen(int **pierce_mat, int arl, int arh, int acl, int ach);
void S_solve(int ** S_ind, double * P_vec, double ** P_mat, int m, int n, double p0, double w, double t); //overwrites vector contents with
void clean_up(int ** S_ind, int ** pierce_mat, double * P_vec, int m, int n, double ** P_mat, double * P_out, double p0, double w, double t_it);
void d_mat_vec_mult_special(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c, int m, int n);
void d_mat_vec_mult_special2(double ** a, int arl, int arh, int acl, int ach, double * b, int brl, int brh, double * c, int m, int n);

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
