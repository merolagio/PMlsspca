#include <RcppEigen.h>
#include <Eigen/Dense>

#include "Cpp_tuto25_LinAlg2.h"

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericVector;

// [[Rcpp::export]]
Eigen::MatrixXd aatC(const Eigen::Map<Eigen::MatrixXd>& A) {
  const int p = A.rows();
  Eigen::MatrixXd M(p, p);
  M.noalias() = A * A.adjoint();
  return M;
}

// [[Rcpp::export]]
Eigen::MatrixXd ataC(const Eigen::Map<Eigen::MatrixXd>& A) {
  const int p = A.cols();
  Eigen::MatrixXd M(p, p);
  M.noalias() = A.adjoint() * A;
  return M;
}

// [[Rcpp::export]]
Eigen::MatrixXd atbC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::MatrixXd>& B) {
  const int p = A.cols();
  const int q = B.cols();
  Eigen::MatrixXd M(p, q);
  M.noalias() = A.adjoint() * B;
  return M;
}

// [[Rcpp::export]]
Eigen::MatrixXd abC(const Eigen::Map<Eigen::MatrixXd>& A,
                    const Eigen::Map<Eigen::MatrixXd>& B) {
  const int p = A.rows();
  const int q = B.cols();
  if (A.cols() != B.rows())
    Rf_error("abC: A and B must be compatible for multiplication A * B");
  Eigen::MatrixXd M(p, q);
  M.noalias() = A * B;
  return M;
}

// [[Rcpp::export]]
Eigen::MatrixXd abtC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::MatrixXd>& B) {
  const int p = A.rows();
  const int q = B.rows();
  if (A.cols() != B.cols())
    Rf_error("abtC: A and B must have the same number of columns");
  Eigen::MatrixXd M(p, q);
  M.noalias() = A * B.adjoint();
  return M;
}

// [[Rcpp::export]]
Eigen::MatrixXd atdaC(const Eigen::Map<Eigen::MatrixXd>& A,
                      const Eigen::Map<Eigen::MatrixXd>& D) {
  const int p = A.cols();
  Eigen::MatrixXd M(p, p);
  M.noalias() = A.adjoint() * D * A;
  return M;
}

// [[Rcpp::export]]
Eigen::MatrixXd atdbC(const Eigen::Map<Eigen::MatrixXd>& A,
                      const Eigen::Map<Eigen::MatrixXd>& D,
                      const Eigen::Map<Eigen::MatrixXd>& B) {
  const int p = A.cols();
  const int q = B.cols();
  Eigen::MatrixXd M(p, q);
  M.noalias() = A.adjoint() * D * B;
  return M;
}

// [[Rcpp::export]]
Eigen::VectorXd avC(const Eigen::Map<Eigen::MatrixXd>& A,
                    const Eigen::Map<Eigen::VectorXd>& v) {
  return A * v;
}

// [[Rcpp::export]]
Eigen::VectorXd vtaC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::VectorXd>& v) {
  // returns (v' A)' as a column vector
  return A.adjoint() * v;
}

// [[Rcpp::export]]
double vtvC(const Eigen::Map<Eigen::VectorXd>& v) {
  return v.squaredNorm();
}

// [[Rcpp::export]]
double vtuC(const Eigen::Map<Eigen::VectorXd>& v,
            const Eigen::Map<Eigen::VectorXd>& u) {
  return v.dot(u);
}

// [[Rcpp::export]]
Eigen::MatrixXd scaleC(const Eigen::Map<Eigen::MatrixXd>& A,
                       bool center = true, bool scale = true) {

  Eigen::MatrixXd C = A; // copy is required (we modify in place)
  Eigen::VectorXd m(A.cols());

  if (center) {
    m = C.colwise().mean();
    C.array().rowwise() -= m.array().transpose();
  }
  if (scale) {
    m = C.colwise().norm();
    C.array().rowwise() /= m.array().transpose();
  }
  return C;
}

// [[Rcpp::export]]
Eigen::MatrixXd scaleColsC(const Eigen::Map<Eigen::MatrixXd>& A, int normtype,
                           const Eigen::Map<Eigen::VectorXd>& sig) {
  // sig is vector of coefficients to multiply each column (usually +/- 1)

  Eigen::MatrixXd M = A; // copy is required (we return a scaled matrix)

  if (normtype == 2) {
    const Eigen::VectorXd l2 = A.colwise().norm();
    M.array().rowwise() /= (l2.array() * sig.array()).transpose();
  } else if (normtype == 1) {
    const Eigen::VectorXd l1 = A.array().abs().colwise().sum().matrix();
    M.array().rowwise() /= (l1.array() * sig.array()).transpose();
  } else {
    M.array().rowwise() /= sig.array().transpose();
  }

  return M;
}

// [[Rcpp::export]]
Rcpp::NumericVector sortC(Rcpp::NumericVector v) {
  Rcpp::NumericVector y = Rcpp::clone(v);
  std::sort(y.begin(), y.end());
  return y;
}

// [[Rcpp::export]]
double sumdiagC(const Eigen::Map<Eigen::MatrixXd>& S) {
  return S.trace();
}


//======================================
//        helpers for lsspca
//======================================

// [[Rcpp::export]]
Rcpp::List EigenC(const Eigen::Map<Eigen::MatrixXd>& S) {
  if ((S - S.adjoint()).norm() > 1e-8 * S.norm())
    Rcpp::stop("S in EigenC must be square");
  if (S != S.adjoint())
    Rcpp::stop("S in EigenC must be symmetric");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  return Rcpp::List::create(
    Rcpp::Named("vec") = Eigen::MatrixXd(es.eigenvectors()).rowwise().reverse(),
    Rcpp::Named("val") = Eigen::VectorXd(es.eigenvalues()).reverse()
  );
}


// [[Rcpp::export]]
Eigen::VectorXd EigenvaluesC(const Eigen::Map<Eigen::MatrixXd>& S) {
  if ((S - S.adjoint()).norm() > 1e-8 * S.norm())
    Rcpp::stop("S in EigenC must be square");
  if (S != S.adjoint())
    Rcpp::stop("S in EigenC must be symmetric");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S, Eigen::EigenvaluesOnly);
  return Eigen::VectorXd(es.eigenvalues().reverse());
}


// [[Rcpp::export]]
Rcpp::List GenEigenC(const Eigen::Map<Eigen::MatrixXd>& A,
                     const Eigen::Map<Eigen::MatrixXd>& B) {
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A, B);
  return Rcpp::List::create(
    Rcpp::Named("vec") = Eigen::MatrixXd(es.eigenvectors()).rowwise().reverse(),
    Rcpp::Named("val") = Eigen::VectorXd(es.eigenvalues()).reverse()
  );
}

// [[Rcpp::export]]
Eigen::MatrixXd solveC(const Eigen::Map<Eigen::MatrixXd>& S) {
  const int d = S.cols();
  Eigen::LLT<Eigen::MatrixXd> llt(S);
  return llt.solve(Eigen::MatrixXd::Identity(d, d));
}
