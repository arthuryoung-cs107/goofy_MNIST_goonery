#ifndef line_search_H  /* Include guard */
#define line_search_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

double line_search_tau(gsl_matrix * A, gsl_matrix * Xk, gsl_vector * bk, gsl_vector * gk, gsl_vector * xvec, gsl_vector * xvec_work, gsl_vector * b_work, double tau_max);

#endif // line_search_H
