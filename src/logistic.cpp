#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <math.h>
#include <stdio.h>

#include "logistic.h"
#include "debug.h"

// I need to bundle all the data that goes to the function to optimze
// together.
typedef struct {
  gsl_matrix_int *X;
  gsl_vector_int *nlev;
  gsl_vector *y;
  gsl_matrix *Xc; // Continuous covariates matrix Nobs x Kc (NULL if not used).
  double lambdaL1;
  double lambdaL2;
} fix_parm_mixed_T;

double fLogit_mixed(gsl_vector *beta, gsl_matrix_int *X, gsl_vector_int *nlev,
                    gsl_matrix *Xc, gsl_vector *y, double lambdaL1,
                    double lambdaL2) {
  int n = y->size;
  int npar = beta->size;
  double total = 0;
  double aux = 0;

  // Changed loop start at 1 instead of 0 to avoid regularization of
  // beta_0*\/
  // #pragma omp parallel for reduction (+:total)
  for (int i = 1; i < npar; ++i)
    total += beta->data[i] * beta->data[i];
  total = (-total * lambdaL2 / 2);
  // #pragma omp parallel for reduction (+:aux)
  for (int i = 1; i < npar; ++i)
    aux += (beta->data[i] > 0 ? beta->data[i] : -beta->data[i]);
  total = total - aux * lambdaL1;
  // #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y)
  // #reduction (+:total)
  for (int i = 0; i < n; ++i) {
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (size_t k = 0; k < X->size2; ++k) {
      if (gsl_matrix_int_get(X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(X, i, k) - 1 + iParm];
      iParm += nlev->data[k] - 1;
    }
    for (size_t k = 0; k < (Xc->size2); ++k)
      Xbetai += gsl_matrix_get(Xc, i, k) * beta->data[iParm++];
    total += y->data[i] * Xbetai - gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
}

void logistic_mixed_pred(gsl_vector *beta,     // Vector of parameters
                                               // length = 1 + Sum_k(C_k -1)
                         gsl_matrix_int *X,    // Matrix Nobs x K
                         gsl_vector_int *nlev, // Vector with number categories
                         gsl_matrix *Xc,       // Continuous covariates matrix:
                                               // obs x Kc (NULL if not used).
                         gsl_vector *yhat) {   // Vector of prob. predicted by
                                               // the logistic
  for (size_t i = 0; i < X->size1; ++i) {
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (size_t k = 0; k < X->size2; ++k) {
      if (gsl_matrix_int_get(X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(X, i, k) - 1 + iParm];
      iParm += nlev->data[k] - 1;
    }
    // Adding the continuous.
    for (size_t k = 0; k < (Xc->size2); ++k)
      Xbetai += gsl_matrix_get(Xc, i, k) * beta->data[iParm++];
    yhat->data[i] = 1 / (1 + gsl_sf_exp(-Xbetai));
  }
}

// The gradient of f, df = (df/dx, df/dy).
void wgsl_mixed_optim_df(const gsl_vector *beta, void *params,
                         gsl_vector *out) {
  fix_parm_mixed_T *p = (fix_parm_mixed_T *)params;
  int n = p->y->size;
  int K = p->X->size2;
  int Kc = p->Xc->size2;
  int npar = beta->size;

  // Intitialize gradient out necessary?
  for (int i = 0; i < npar; ++i)
    out->data[i] = 0;

  // Changed loop start at 1 instead of 0 to avoid regularization of beta 0.
  for (int i = 1; i < npar; ++i)
    out->data[i] = p->lambdaL2 * beta->data[i];
  for (int i = 1; i < npar; ++i)
    out->data[i] += p->lambdaL1 * ((beta->data[i] > 0) - (beta->data[i] < 0));

  for (int i = 0; i < n; ++i) {
    double pn = 0;
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm];
      iParm += p->nlev->data[k] - 1;
    }

    // Adding the continuous.
    for (int k = 0; k < Kc; ++k)
      Xbetai += gsl_matrix_get(p->Xc, i, k) * beta->data[iParm++];

    pn = -(p->y->data[i] - 1 / (1 + gsl_sf_exp(-Xbetai)));

    out->data[0] += pn;
    iParm = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        out->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm] += pn;
      iParm += p->nlev->data[k] - 1;
    }

    // Adding the continuous.
    for (int k = 0; k < Kc; ++k) {
      out->data[iParm++] += gsl_matrix_get(p->Xc, i, k) * pn;
    }
  }
}

// The Hessian of f.
void wgsl_mixed_optim_hessian(const gsl_vector *beta, void *params,
                              gsl_matrix *out) {
  fix_parm_mixed_T *p = (fix_parm_mixed_T *)params;
  int n = p->y->size;
  int K = p->X->size2;
  int Kc = p->Xc->size2;
  int npar = beta->size;
  gsl_vector *gn = gsl_vector_safe_alloc(npar); // gn

  // Intitialize Hessian out necessary ???
  gsl_matrix_set_zero(out);

  /* Changed loop start at 1 instead of 0 to avoid regularization of beta 0*/
  for (int i = 1; i < npar; ++i)
    gsl_matrix_set(out, i, i, (p->lambdaL2)); // Double-check this.

  // L1 penalty not working yet, as not differentiable, I may need to
  // do coordinate descent (as in glm_net)
  for (int i = 0; i < n; ++i) {
    double pn = 0;
    double aux = 0;
    double Xbetai = beta->data[0];
    int iParm1 = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm1];
      iParm1 += p->nlev->data[k] - 1; //-1?
    }

    // Adding the continuous.
    for (int k = 0; k < Kc; ++k)
      Xbetai += gsl_matrix_get(p->Xc, i, k) * beta->data[iParm1++];

    pn = 1 / (1 + gsl_sf_exp(-Xbetai));

    // Add a protection for pn very close to 0 or 1?
    aux = pn * (1 - pn);

    // Calculate sub-gradient vector gn.
    gsl_vector_set_zero(gn);
    gn->data[0] = 1;
    iParm1 = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        gn->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm1] = 1;
      iParm1 += p->nlev->data[k] - 1;
    }

    // Adding the continuous.
    for (int k = 0; k < Kc; ++k) {
      gn->data[iParm1++] = gsl_matrix_get(p->Xc, i, k);
    }

    for (int k1 = 0; k1 < npar; ++k1)
      if (gn->data[k1] != 0)
        for (int k2 = 0; k2 < npar; ++k2)
          if (gn->data[k2] != 0)
            *gsl_matrix_ptr(out, k1, k2) += (aux * gn->data[k1] * gn->data[k2]);
  }
  gsl_vector_free(gn);
}

double wgsl_mixed_optim_f(gsl_vector *v, void *params) {
  fix_parm_mixed_T *p = (fix_parm_mixed_T *)params;
  return fLogit_mixed(v, p->X, p->nlev, p->Xc, p->y, p->lambdaL1, p->lambdaL2);
}

// Compute both f and df together.
void wgsl_mixed_optim_fdf(gsl_vector *x, void *params, double *f,
                          gsl_vector *df) {
  *f = wgsl_mixed_optim_f(x, params);
  wgsl_mixed_optim_df(x, params, df);
}

// Xc is the matrix of continuous covariates, Nobs x Kc (NULL if not used).
int logistic_mixed_fit(gsl_vector *beta, gsl_matrix_int *X,
                       gsl_vector_int *nlev, gsl_matrix *Xc, gsl_vector *y,
                       double lambdaL1, double lambdaL2) {
  // double mLogLik = 0;
  fix_parm_mixed_T p;
  int npar = beta->size;
  int iter = 0;
  double maxchange = 0;

  // Intializing fix parameters.
  p.X = X;
  p.Xc = Xc;
  p.nlev = nlev;
  p.y = y;
  p.lambdaL1 = lambdaL1;
  p.lambdaL2 = lambdaL2;

  // Initial fit.
  // auto mLogLik = wgsl_mixed_optim_f(beta, &p);

  gsl_matrix *myH = gsl_matrix_safe_alloc(npar, npar); // Hessian matrix.
  gsl_vector *stBeta = gsl_vector_safe_alloc(npar);    // Direction to move.

  gsl_vector *myG = gsl_vector_safe_alloc(npar); // Gradient.
  gsl_vector *tau = gsl_vector_safe_alloc(npar); // tau for QR.

  for (iter = 0; iter < 100; iter++) {
    wgsl_mixed_optim_hessian(beta, &p, myH); // Calculate Hessian.
    wgsl_mixed_optim_df(beta, &p, myG);      // Calculate Gradient.
    gsl_linalg_QR_decomp(myH, tau);          // Calculate next beta.
    gsl_linalg_QR_solve(myH, tau, myG, stBeta);
    gsl_vector_sub(beta, stBeta);

    // Monitor convergence.
    maxchange = 0;
    for (int i = 0; i < npar; i++)
      if (maxchange < fabs(stBeta->data[i]))
        maxchange = fabs(stBeta->data[i]);

    if (maxchange < 1E-4)
      break;
  }

  // Final fit.
  // mLogLik = wgsl_mixed_optim_f(beta, &p);

  gsl_vector_free(tau);
  gsl_vector_free(stBeta);
  gsl_vector_free(myG);
  gsl_matrix_free(myH);

  return 0;
}

/***************/
/* Categorical */
/***************/

// I need to bundle all the data that goes to the function to optimze
// together.
typedef struct {
  gsl_matrix_int *X;
  gsl_vector_int *nlev;
  gsl_vector *y;
  double lambdaL1;
  double lambdaL2;
} fix_parm_cat_T;

double fLogit_cat(gsl_vector *beta, gsl_matrix_int *X, gsl_vector_int *nlev,
                  gsl_vector *y, double lambdaL1, double lambdaL2) {
  int n = y->size;
  int npar = beta->size;
  double total = 0;
  double aux = 0;

  // omp_set_num_threads(ompthr); /\* Changed loop start at 1 instead
  // of 0 to avoid regularization of beta 0*\/ /\*#pragma omp parallel
  // for reduction (+:total)*\/
  for (int i = 1; i < npar; ++i)
    total += beta->data[i] * beta->data[i];
  total = (-total * lambdaL2 / 2);

  // /\*#pragma omp parallel for reduction (+:aux)*\/
  for (int i = 1; i < npar; ++i)
    aux += (beta->data[i] > 0 ? beta->data[i] : -beta->data[i]);
  total = total - aux * lambdaL1;

  // #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y)
  // #reduction (+:total)
  for (int i = 0; i < n; ++i) {
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (size_t k = 0; k < X->size2; ++k) {
      if (gsl_matrix_int_get(X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(X, i, k) - 1 + iParm];
      iParm += nlev->data[k] - 1;
    }
    total += y->data[i] * Xbetai - gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
}

void logistic_cat_pred(gsl_vector *beta,     // Vector of parameters
                                             // length = 1 + Sum_k(C_k-1).
                       gsl_matrix_int *X,    // Matrix Nobs x K
                       gsl_vector_int *nlev, // Vector with #categories
                       gsl_vector *yhat) {   // Vector of prob. predicted by
                                             // the logistic.
  for (size_t i = 0; i < X->size1; ++i) {
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (size_t k = 0; k < X->size2; ++k) {
      if (gsl_matrix_int_get(X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(X, i, k) - 1 + iParm];
      iParm += nlev->data[k] - 1;
    }
    yhat->data[i] = 1 / (1 + gsl_sf_exp(-Xbetai));
  }
}

// The gradient of f, df = (df/dx, df/dy).
void wgsl_cat_optim_df(const gsl_vector *beta, void *params, gsl_vector *out) {
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  int n = p->y->size;
  int K = p->X->size2;
  int npar = beta->size;

  // Intitialize gradient out necessary?
  for (int i = 0; i < npar; ++i)
    out->data[i] = 0;

  // Changed loop start at 1 instead of 0 to avoid regularization of beta 0.
  for (int i = 1; i < npar; ++i)
    out->data[i] = p->lambdaL2 * beta->data[i];
  for (int i = 1; i < npar; ++i)
    out->data[i] += p->lambdaL1 * ((beta->data[i] > 0) - (beta->data[i] < 0));

  for (int i = 0; i < n; ++i) {
    double pn = 0;
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm];
      iParm += p->nlev->data[k] - 1;
    }

    pn = -(p->y->data[i] - 1 / (1 + gsl_sf_exp(-Xbetai)));

    out->data[0] += pn;
    iParm = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        out->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm] += pn;
      iParm += p->nlev->data[k] - 1;
    }
  }
}

// The Hessian of f.
void wgsl_cat_optim_hessian(const gsl_vector *beta, void *params,
                            gsl_matrix *out) {
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  int n = p->y->size;
  int K = p->X->size2;
  int npar = beta->size;

  // Intitialize Hessian out necessary.
  gsl_matrix_set_zero(out);

  // Changed loop start at 1 instead of 0 to avoid regularization of beta.
  for (int i = 1; i < npar; ++i)
    gsl_matrix_set(out, i, i, (p->lambdaL2)); // Double-check this.

  // L1 penalty not working yet, as not differentiable, I may need to
  // do coordinate descent (as in glm_net).
  for (int i = 0; i < n; ++i) {
    double pn = 0;
    double aux = 0;
    double Xbetai = beta->data[0];
    int iParm2 = 1;
    int iParm1 = 1;
    for (int k = 0; k < K; ++k) {
      if (gsl_matrix_int_get(p->X, i, k) > 0)
        Xbetai += beta->data[gsl_matrix_int_get(p->X, i, k) - 1 + iParm1];
      iParm1 += p->nlev->data[k] - 1; //-1?
    }

    pn = 1 / (1 + gsl_sf_exp(-Xbetai));

    // Add a protection for pn very close to 0 or 1?
    aux = pn * (1 - pn);
    *gsl_matrix_ptr(out, 0, 0) += aux;
    iParm2 = 1;
    for (int k2 = 0; k2 < K; ++k2) {
      if (gsl_matrix_int_get(p->X, i, k2) > 0)
        *gsl_matrix_ptr(out, 0, gsl_matrix_int_get(p->X, i, k2) - 1 + iParm2) +=
            aux;
      iParm2 += p->nlev->data[k2] - 1; //-1?
    }
    iParm1 = 1;
    for (int k1 = 0; k1 < K; ++k1) {
      if (gsl_matrix_int_get(p->X, i, k1) > 0)
        *gsl_matrix_ptr(out, gsl_matrix_int_get(p->X, i, k1) - 1 + iParm1, 0) +=
            aux;
      iParm2 = 1;
      for (int k2 = 0; k2 < K; ++k2) {
        if ((gsl_matrix_int_get(p->X, i, k1) > 0) &&
            (gsl_matrix_int_get(p->X, i, k2) > 0))
          *gsl_matrix_ptr(out, gsl_matrix_int_get(p->X, i, k1) - 1 + iParm1,
                          gsl_matrix_int_get(p->X, i, k2) - 1 + iParm2) += aux;
        iParm2 += p->nlev->data[k2] - 1; //-1?
      }
      iParm1 += p->nlev->data[k1] - 1; //-1?
    }
  }
}

double wgsl_cat_optim_f(gsl_vector *v, void *params) {
  double mLogLik = 0;
  fix_parm_cat_T *p = (fix_parm_cat_T *)params;
  mLogLik = fLogit_cat(v, p->X, p->nlev, p->y, p->lambdaL1, p->lambdaL2);
  return mLogLik;
}

// Compute both f and df together.
void wgsl_cat_optim_fdf(gsl_vector *x, void *params, double *f,
                        gsl_vector *df) {
  *f = wgsl_cat_optim_f(x, params);
  wgsl_cat_optim_df(x, params, df);
}

int logistic_cat_fit(gsl_vector *beta, gsl_matrix_int *X, gsl_vector_int *nlev,
                     gsl_vector *y, double lambdaL1, double lambdaL2) {
  // double mLogLik = 0;
  fix_parm_cat_T p;
  int npar = beta->size;
  int iter = 0;
  double maxchange = 0;

  // Intializing fix parameters.
  p.X = X;
  p.nlev = nlev;
  p.y = y;
  p.lambdaL1 = lambdaL1;
  p.lambdaL2 = lambdaL2;

#ifdef _RPR_DEBUG_
  // Initial fit.
  auto mLogLik = wgsl_cat_optim_f(beta, &p);
#endif

  gsl_matrix *myH = gsl_matrix_safe_alloc(npar, npar); // Hessian matrix.
  gsl_vector *stBeta = gsl_vector_safe_alloc(npar);    // Direction to move.

  gsl_vector *myG = gsl_vector_safe_alloc(npar); // Gradient.
  gsl_vector *tau = gsl_vector_safe_alloc(npar); // tau for QR.

  for (iter = 0; iter < 100; iter++) {
    wgsl_cat_optim_hessian(beta, &p, myH); // Calculate Hessian.
    wgsl_cat_optim_df(beta, &p, myG);      // Calculate Gradient.
    gsl_linalg_QR_decomp(myH, tau);        // Calculate next beta.
    gsl_linalg_QR_solve(myH, tau, myG, stBeta);
    gsl_vector_sub(beta, stBeta);

    // Monitor convergence.
    maxchange = 0;
    for (int i = 0; i < npar; i++)
      if (maxchange < fabs(stBeta->data[i]))
        maxchange = fabs(stBeta->data[i]);

#ifdef _RPR_DEBUG_
    mLogLik = wgsl_cat_optim_f(beta, &p);
#endif

    if (maxchange < 1E-4)
      break;
  }

  // Final fit.
  // mLogLik = wgsl_cat_optim_f(beta, &p);

  gsl_vector_free(tau);
  gsl_vector_free(stBeta);
  gsl_vector_free(myG);
  gsl_matrix_free(myH);

  return 0;
}

/***************/
/* Continuous  */
/***************/

// I need to bundle all the data that goes to the function to optimze
// together.
typedef struct {
  gsl_matrix *Xc; // continuous covariates; Matrix Nobs x Kc
  gsl_vector *y;
  double lambdaL1;
  double lambdaL2;
} fix_parm_cont_T;

double fLogit_cont(const gsl_vector *beta, const gsl_matrix *Xc, const gsl_vector *y,
                   double lambdaL1, double lambdaL2) {
  int n = y->size;
  int npar = beta->size;
  double total = 0;
  double aux = 0;

  // omp_set_num_threads(ompthr); /\* Changed loop start at 1 instead
  // of 0 to avoid regularization of beta_0*\/ /\*#pragma omp parallel
  // for reduction (+:total)*\/
  for (int i = 1; i < npar; ++i)
    total += beta->data[i] * beta->data[i];
  total = (-total * lambdaL2 / 2);

  // /\*#pragma omp parallel for reduction (+:aux)*\/
  for (int i = 1; i < npar; ++i)
    aux += (beta->data[i] > 0 ? beta->data[i] : -beta->data[i]);
  total = total - aux * lambdaL1;

  // #pragma omp parallel for schedule(static) shared(n,beta,X,nlev,y)
  // #reduction (+:total)
  for (int i = 0; i < n; ++i) {
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (size_t k = 0; k < (Xc->size2); ++k)
      Xbetai += gsl_matrix_get(Xc, i, k) * beta->data[iParm++];
    total += y->data[i] * Xbetai - gsl_sf_log_1plusx(gsl_sf_exp(Xbetai));
  }
  return -total;
}

void logistic_cont_pred(gsl_vector *beta,   // Vector of parameters
                                            // length = 1 + Sum_k(C_k-1).
                        gsl_matrix *Xc,     // Continuous covariates matrix,
                                            // Nobs x Kc (NULL if not used).
                        gsl_vector *yhat) { // Vector of prob. predicted by
                                            // the logistic.
  for (size_t i = 0; i < Xc->size1; ++i) {
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (size_t k = 0; k < (Xc->size2); ++k)
      Xbetai += gsl_matrix_get(Xc, i, k) * beta->data[iParm++];
    yhat->data[i] = 1 / (1 + gsl_sf_exp(-Xbetai));
  }
}

// The gradient of f, df = (df/dx, df/dy).
void wgsl_cont_optim_df(const gsl_vector *beta, const void *params, gsl_vector *out) {
  fix_parm_cont_T *p = (fix_parm_cont_T *)params;
  int n = p->y->size;
  int Kc = p->Xc->size2;
  int npar = beta->size;

  // Intitialize gradient out necessary?
  for (int i = 0; i < npar; ++i)
    out->data[i] = 0;

  // Changed loop start at 1 instead of 0 to avoid regularization of beta 0.
  for (int i = 1; i < npar; ++i)
    out->data[i] = p->lambdaL2 * beta->data[i];
  for (int i = 1; i < npar; ++i)
    out->data[i] += p->lambdaL1 * ((beta->data[i] > 0) - (beta->data[i] < 0));

  for (int i = 0; i < n; ++i) {
    double pn = 0;
    double Xbetai = beta->data[0];
    int iParm = 1;
    for (int k = 0; k < Kc; ++k)
      Xbetai += gsl_matrix_get(p->Xc, i, k) * beta->data[iParm++];

    pn = -(p->y->data[i] - 1 / (1 + gsl_sf_exp(-Xbetai)));

    out->data[0] += pn;
    iParm = 1;

    // Adding the continuous.
    for (int k = 0; k < Kc; ++k) {
      out->data[iParm++] += gsl_matrix_get(p->Xc, i, k) * pn;
    }
  }
}

// The Hessian of f.
void wgsl_cont_optim_hessian(const gsl_vector *beta, void *params,
                             gsl_matrix *out) {
  fix_parm_cont_T *p = (fix_parm_cont_T *)params;
  int n = p->y->size;
  int Kc = p->Xc->size2;
  int npar = beta->size;
  gsl_vector *gn = gsl_vector_safe_alloc(npar); // gn.

  // Intitialize Hessian out necessary ???

  gsl_matrix_set_zero(out);

  // Changed loop start at 1 instead of 0 to avoid regularization of
  // beta 0.
  for (int i = 1; i < npar; ++i)
    gsl_matrix_set(out, i, i, (p->lambdaL2)); // Double-check this.

  // L1 penalty not working yet, as not differentiable, I may need to
  // do coordinate descent (as in glm_net).
  for (int i = 0; i < n; ++i) {
    double pn = 0;
    double aux = 0;
    double Xbetai = beta->data[0];
    int iParm1 = 1;
    for (int k = 0; k < Kc; ++k)
      Xbetai += gsl_matrix_get(p->Xc, i, k) * beta->data[iParm1++];

    pn = 1 / (1 + gsl_sf_exp(-Xbetai));

    // Add a protection for pn very close to 0 or 1?
    aux = pn * (1 - pn);

    // Calculate sub-gradient vector gn.
    gsl_vector_set_zero(gn);
    gn->data[0] = 1;
    iParm1 = 1;
    for (int k = 0; k < Kc; ++k) {
      gn->data[iParm1++] = gsl_matrix_get(p->Xc, i, k);
    }

    for (int k1 = 0; k1 < npar; ++k1)
      if (gn->data[k1] != 0)
        for (int k2 = 0; k2 < npar; ++k2)
          if (gn->data[k2] != 0)
            *gsl_matrix_ptr(out, k1, k2) += (aux * gn->data[k1] * gn->data[k2]);
  }
  gsl_vector_free(gn);
}

double wgsl_cont_optim_f(const gsl_vector *v, const void *params) {
  double mLogLik = 0;
  fix_parm_cont_T *p = (fix_parm_cont_T *)params;
  mLogLik = fLogit_cont(v, p->Xc, p->y, p->lambdaL1, p->lambdaL2);
  return mLogLik;
}

// Compute both f and df together.
void wgsl_cont_optim_fdf(const gsl_vector *x, const void *params, double *f,
                         gsl_vector *df) {
  *f = wgsl_cont_optim_f(x, params);
  wgsl_cont_optim_df(x, params, df);
}

int logistic_cont_fit(gsl_vector *beta,
                      gsl_matrix *Xc, // Continuous covariates matrix,
                                      // Nobs x Kc (NULL if not used).
                      gsl_vector *y, double lambdaL1, double lambdaL2) {

  fix_parm_cont_T p;
  int npar = beta->size;
  int iter = 0;
  double maxchange = 0;

  // Initializing fix parameters.
  p.Xc = Xc;
  p.y = y;
  p.lambdaL1 = lambdaL1;
  p.lambdaL2 = lambdaL2;

#ifdef _RPR_DEBUG_
  // Initial fit.
  auto mLogLik = wgsl_cont_optim_f(beta, &p);
#endif

  gsl_matrix *myH = gsl_matrix_safe_alloc(npar, npar); // Hessian matrix.
  gsl_vector *stBeta = gsl_vector_safe_alloc(npar);    // Direction to move.

  gsl_vector *myG = gsl_vector_safe_alloc(npar); // Gradient.
  gsl_vector *tau = gsl_vector_safe_alloc(npar); // tau for QR.

  for (iter = 0; iter < 100; iter++) {
    wgsl_cont_optim_hessian(beta, &p, myH); // Calculate Hessian.
    wgsl_cont_optim_df(beta, &p, myG);      // Calculate Gradient.
    gsl_linalg_QR_decomp(myH, tau);         // Calculate next beta.
    gsl_linalg_QR_solve(myH, tau, myG, stBeta);
    gsl_vector_sub(beta, stBeta);

    // Monitor convergence.
    maxchange = 0;
    for (int i = 0; i < npar; i++)
      if (maxchange < fabs(stBeta->data[i]))
        maxchange = fabs(stBeta->data[i]);

#ifdef _RPR_DEBUG_
    mLogLik = wgsl_cont_optim_f(beta, &p);
#endif

    if (maxchange < 1E-4)
      break;
  }

  // Final fit.
  // mLogLik = wgsl_cont_optim_f(beta, &p);

  gsl_vector_free(tau);
  gsl_vector_free(stBeta);
  gsl_vector_free(myG);
  gsl_matrix_free(myH);

  return 0;
}
