#ifndef GSL_tomfoolery_H  /* Include guard */
#define GSL_tomfoolery_H

#include "GSL_tomfoolery.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

void NRmat_2_GSLmat(double ** matrix, int arl, int arh, int acl, int ach, gsl_matrix * GSL_mat);
void GSLmat_2_NRmat(double ** matrix, int arl, int arh, int acl, int ach, gsl_matrix * GSL_mat);
void GSLvec_2_NRvec(double * vec, int arl, int arh, gsl_vector * GSL_vec);
void NRvec_2_GSLvec(double * vec, int arl, int arh, gsl_vector * GSL_vec);

void GSL_mat2vec(gsl_matrix * mat, gsl_vector * vec);
void GSL_vec2mat(gsl_matrix * mat, gsl_vector * vec);
void G_func_GSL(gsl_matrix * A, gsl_matrix * X, gsl_vector * b, gsl_vector * b_work, gsl_vector * x_vec);
void matrix_shrink_GSL(gsl_matrix * X, gsl_matrix * U, gsl_vector * S, gsl_matrix * V, gsl_vector * work, gsl_matrix * S_work, gsl_matrix * nn_work, double nu);
void vector_shrink_GSL(gsl_vector * S, double nu);


#endif // GSL_tomfoolery_H
