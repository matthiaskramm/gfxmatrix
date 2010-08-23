/* gsltools.c 

   GSL helper functions.

   Copyright (c) 2007,2008,2009,2010 Matthias Kramm <kramm@quiss.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include "gsltools.h"

int my_gsl_linalg_LU_decomp (gsl_matrix * A, gsl_permutation * p, int *signum)
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("LU decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != A->size1)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      size_t i, j, k;

      *signum = 1;
      gsl_permutation_init (p);

      for (j = 0; j < N - 1; j++)
        {
          /* Find maximum in the j-th column */

          double ajj, max = fabs (gsl_matrix_get (A, j, j));
          size_t i_pivot = j;

          for (i = j + 1; i < N; i++)
            {
              double aij = fabs (gsl_matrix_get (A, i, j));

              if (aij > max)
                {
                  max = aij;
                  i_pivot = i;
                }
            }

          if (i_pivot != j)
            {
              gsl_matrix_swap_rows (A, j, i_pivot);
              gsl_permutation_swap (p, j, i_pivot);
              *signum = -(*signum);
            }

          ajj = gsl_matrix_get (A, j, j);

          if (fabs(ajj) > 0.0000001)
            {
              //printf("Using %f (column %d, row %d) as pivot\n", ajj, j, i_pivot);

              for (i = j + 1; i < N; i++)
                {
                  double aij = gsl_matrix_get (A, i, j) / ajj;
                  gsl_matrix_set (A, i, j, aij);

                  for (k = j + 1; k < N; k++)
                    {
                      double aik = gsl_matrix_get (A, i, k);
                      double ajk = gsl_matrix_get (A, j, k);
                      //printf("Setting element %d,%d to %f\n", i,k, aik - aij * ajk);
                      gsl_matrix_set (A, i, k, aik - aij * ajk);
                    }
                }
            } 
            else 
            {
                fprintf(stderr, "GSL Error: Pivot is zero (column %d, row %d/%d)\n", j, i_pivot, p->data[i_pivot]);
                return -1;
            }
        }
      
      return GSL_SUCCESS;
    }
}

/* Handle the general case of the SVD with tolerance and rank */

int
my_gsl_multifit_linear_svd (const gsl_matrix * X,
                         const gsl_vector * y,
                         double tol,
                         size_t * rank,
                         gsl_vector * c,
                         gsl_matrix * cov,
                         double *chisq, gsl_multifit_linear_workspace * work)
{
  if (X->size1 != y->size)
    {
      GSL_ERROR
        ("number of observations in y does not match rows of matrix X",
         GSL_EBADLEN);
    }
  else if (X->size2 != c->size)
    {
      GSL_ERROR ("number of parameters c does not match columns of matrix X",
                 GSL_EBADLEN);
    }
  else if (cov->size1 != cov->size2)
    {
      GSL_ERROR ("covariance matrix is not square", GSL_ENOTSQR);
    }
  else if (c->size != cov->size1)
    {
      GSL_ERROR
        ("number of parameters does not match size of covariance matrix",
         GSL_EBADLEN);
    }
  else if (X->size1 != work->n || X->size2 != work->p)
    {
      GSL_ERROR
        ("size of workspace does not match size of observation matrix",
         GSL_EBADLEN);
    }
  else if (tol <= 0)
    {
      GSL_ERROR ("tolerance must be positive", GSL_EINVAL);
    }
  else
    {
      const size_t n = X->size1;
      const size_t p = X->size2;

      size_t i, j, p_eff;

      gsl_matrix *A = work->A;
      gsl_matrix *Q = work->Q;
      gsl_matrix *QSI = work->QSI;
      gsl_vector *S = work->S;
      gsl_vector *xt = work->xt;
      gsl_vector *D = work->D;

      /* Copy X to workspace,  A <= X */

      gsl_matrix_memcpy (A, X);

      /* Balance the columns of the matrix A */

      gsl_linalg_balance_columns (A, D);

      /* Decompose A into U S Q^T */
      gsl_linalg_SV_decomp_mod (A, QSI, Q, S, xt);

      /* Solve y = A c for c */

      gsl_blas_dgemv (CblasTrans, 1.0, A, y, 0.0, xt);

      /* Scale the matrix Q,  Q' = Q S^-1 */
      gsl_matrix_memcpy (QSI, Q);

      {
        double alpha0 = gsl_vector_get (S, 0);
        p_eff = 0;

        for (j = 0; j < p; j++)
          {
            gsl_vector_view column = gsl_matrix_column (QSI, j);
            double alpha = gsl_vector_get (S, j);

            if (alpha <= tol * alpha0) {
              alpha = 0.0;
            } else {
              alpha = 1.0 / alpha;
              p_eff++;
            }

            gsl_vector_scale (&column.vector, alpha);
          }

        *rank = p_eff;
      }

      gsl_vector_set_zero (c);

      gsl_blas_dgemv (CblasNoTrans, 1.0, QSI, xt, 0.0, c);

      /* Unscale the balancing factors */

      gsl_vector_div (c, D);

      /* Compute chisq, from residual r = y - X c */

      {
        double s2 = 0, r2 = 0;

        for (i = 0; i < n; i++)
          {
            double yi = gsl_vector_get (y, i);
            gsl_vector_const_view row = gsl_matrix_const_row (X, i);
            double y_est, ri;
            gsl_blas_ddot (&row.vector, c, &y_est);
            ri = yi - y_est;
            r2 += ri * ri;
          }

        s2 = r2 / (n - p_eff);   /* p_eff == rank */

        *chisq = r2;

        /* Form variance-covariance matrix cov = s2 * (Q S^-1) (Q S^-1)^T */

        for (i = 0; i < p; i++)
          {
            gsl_vector_view row_i = gsl_matrix_row (QSI, i);
            double d_i = gsl_vector_get (D, i);

            for (j = i; j < p; j++)
              {
                gsl_vector_view row_j = gsl_matrix_row (QSI, j);
                double d_j = gsl_vector_get (D, j);
                double s;

                gsl_blas_ddot (&row_i.vector, &row_j.vector, &s);

                gsl_matrix_set (cov, i, j, s * s2 / (d_i * d_j));
                gsl_matrix_set (cov, j, i, s * s2 / (d_i * d_j));
              }
          }
      }

      return GSL_SUCCESS;
    }
}

int
my_gsl_multifit_linear (const gsl_matrix * X,
                     const gsl_vector * y,
                     gsl_vector * c,
                     gsl_matrix * cov,
                     double *chisq, gsl_multifit_linear_workspace * work)
{
  size_t rank;
  int status  = my_gsl_multifit_linear_svd (X, y, GSL_DBL_EPSILON, &rank, c,
                                         cov, chisq, work);
  return status;
}

