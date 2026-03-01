#ifndef CPP_TUTO25_LINALG2_H
#define CPP_TUTO25_LINALG2_H

#include <RcppEigen.h>
#include <Eigen/Dense>

// Basic linear algebra helpers (no-copy inputs via Eigen::Map)

// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd aatC(const Eigen::Map<Eigen::MatrixXd>& A);
Eigen::MatrixXd ataC(const Eigen::Map<Eigen::MatrixXd>& A);
Eigen::MatrixXd atbC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::MatrixXd>& B);
Eigen::MatrixXd abC(const Eigen::Map<Eigen::MatrixXd>& A,
                    const Eigen::Map<Eigen::MatrixXd>& B);
Eigen::MatrixXd abtC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::MatrixXd>& B);
Eigen::MatrixXd atdaC(const Eigen::Map<Eigen::MatrixXd>& A,
                      const Eigen::Map<Eigen::MatrixXd>& D);
Eigen::MatrixXd atdbC(const Eigen::Map<Eigen::MatrixXd>& A,
                      const Eigen::Map<Eigen::MatrixXd>& D,
                      const Eigen::Map<Eigen::MatrixXd>& B);

Eigen::VectorXd avC(const Eigen::Map<Eigen::MatrixXd>& A,
                    const Eigen::Map<Eigen::VectorXd>& v);
Eigen::VectorXd vtaC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::VectorXd>& v);

double vtvC(const Eigen::Map<Eigen::VectorXd>& v);
double vtuC(const Eigen::Map<Eigen::VectorXd>& v,
            const Eigen::Map<Eigen::VectorXd>& u);

Eigen::MatrixXd scaleC(const Eigen::Map<Eigen::MatrixXd>& A,
                       bool center, bool scale);

Eigen::MatrixXd scaleColsC(const Eigen::Map<Eigen::MatrixXd>& A, int normtype,
                           const Eigen::Map<Eigen::VectorXd>& sig);

Rcpp::NumericVector sortC(Rcpp::NumericVector v);

double sumdiagC(const Eigen::Map<Eigen::MatrixXd>& S);

// Helpers for LSSPCA
Rcpp::List EigenC(const Eigen::Map<Eigen::MatrixXd>& S);
Rcpp::List GenEigenC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::MatrixXd>& B);
Eigen::MatrixXd solveC(const Eigen::Map<Eigen::MatrixXd>& S);

#endif // CPP_TUTO25_LINALG2_H
