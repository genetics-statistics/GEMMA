#ifndef LOGISTIC_H_
#define LOGISTIC_H_

// Mixed interface.
void logistic_mixed_pred(gsl_vector *beta,     // Vector of parameters
                                               // length = 1+Sum_k(C_k-1)+Kc.
                         gsl_matrix_int *X,    // Matrix Nobs x K.
                         gsl_vector_int *nlev, // Vector with num. categories.
                         gsl_matrix *Xc,       // Continuous covariates matrix
                                               // Nobs x Kc
                         gsl_vector *yhat);    // Vector of prob. predicted by
                                               // the logistic.

int logistic_mixed_fit(gsl_vector *beta,     // Vector of parameters
                                             // length = 1+Sum_k(C_k-1)+Kc
                       gsl_matrix_int *X,    // Matrix Nobs x K.
                       gsl_vector_int *nlev, // Vector with number categories.
                       gsl_matrix *Xc,       // Continuous covariates
                                             // matrix Nobs x Kc
                       gsl_vector *y,        // Vector of prob. to predict.
                       double lambdaL1,      // Reg. L1 0.0 if not used.
                       double lambdaL2);     // Reg. L2 0.0 if not used.

double fLogit_mixed(gsl_vector *beta, gsl_matrix_int *X, gsl_vector_int *nlev,
                    gsl_matrix *Xc, // continuous covariates matrix Nobs x Kc
                    gsl_vector *y, double lambdaL1, double lambdaL2);

// Categorical-only interface.
void logistic_cat_pred(gsl_vector *beta,     // Vector of parameters
                                             // length = 1+Sum_k(C_k-1)+Kc.
                       gsl_matrix_int *X,    // Matrix Nobs x K.
                       gsl_vector_int *nlev, // Vector with number categories.
                       gsl_vector *yhat);    // Vector of prob. predicted by
                                             // the logistic.

int logistic_cat_fit(gsl_vector *beta,     // Vector of parameters
                                           // length = 1+Sum_k(C_k-1)+Kc.
                     gsl_matrix_int *X,    // Matrix Nobs x K .
                     gsl_vector_int *nlev, // Vector with number categories.
                     gsl_vector *y,        // Vector of prob. to predict.
                     double lambdaL1,      // Regularization L1, 0 if not used
                     double lambdaL2);     // Regularization L2, 0 if not used

double fLogit_cat(gsl_vector *beta, gsl_matrix_int *X, gsl_vector_int *nlev,
                  gsl_vector *y, double lambdaL1, double lambdaL2);

// Continuous-only interface.
void logistic_cont_pred(gsl_vector *beta,  // Vector of parameters
                                           // length = 1 + Sum_k(C_k-1) + Kc.
                        gsl_matrix *Xc,    // Continuous cov's matrix Nobs x Kc.
                        gsl_vector *yhat); // Vector of prob. predicted
                                           // by the logistic.

int logistic_cont_fit(gsl_vector *beta, // Vector of parameters
                                        // length = 1+Sum_k(C_k-1)+Kc.
                      gsl_matrix *Xc,   // Continuous cov's matrix Nobs x Kc.
                      gsl_vector *y,    // Vector of prob. to predict.
                      double lambdaL1,  // Regularization L1, 0 if not used.
                      double lambdaL2); // Regularization L2, 0 if not used.

double fLogit_cont(const gsl_vector *beta,
                   const gsl_matrix *Xc, // Continuous covariates matrix Nobs x Kc.
                   const gsl_vector *y, double lambdaL1, double lambdaL2);

#endif
