#include <cmath>
#include <RcppEigen.h>
#include <Eigen/Dense>

// Updated linear-algebra helpers (no-copy interfaces)
#include "Cpp_tuto25_LinAlg2.h"

// [[Rcpp::depends(RcppEigen)]]

using Rcpp::List;
using Rcpp::Named;

//====================================================
// Deflation and explained-variance helpers
//====================================================

// [[Rcpp::export]]
List deflXCforR(const Eigen::Map<Eigen::VectorXd>& comp,
               const Eigen::Map<Eigen::MatrixXd>& K) {
  // Deflates data matrix K by the (possibly non-unit-length) component comp.
  // Returns a list with:
  //   K    = K - (t t' / (t' t)) K
  //   vexp = || (t t' / (t' t)) K ||_F^2

  if (K.rows() != comp.size()) {
    Rf_error("deflXCforR: matrix and component are not compatible for deflation");
  }

  const double tt = comp.squaredNorm();
  if (!std::isfinite(tt) || tt <= 0.0) {
    Rf_error("deflXCforR: component has zero (or non-finite) norm");
  }

  // Avoid forming T = (t t') / (t' t) (dense n x n). Use: T K = t (t' K) / (t' t)
  const Eigen::RowVectorXd ck = comp.transpose() * K;         // 1 x p
  Eigen::MatrixXd M = (comp * ck) / tt;                       // n x p

  const double vexp = M.squaredNorm();
  Eigen::MatrixXd Kdefl = K - M;

  return List::create(Named("K") = Kdefl, Named("vexp") = vexp);
}

// [[Rcpp::export]]
List makeVexpC(const Eigen::Map<Eigen::MatrixXd>& A,
              const Eigen::Map<Eigen::MatrixXd>& S) {
  // Computes marginal and cumulative explained variance for components A
  // relative to covariance/correlation matrix S.

  const int d = A.cols();
  Eigen::VectorXd vexp(d);
  Eigen::VectorXd cvexp(d);

  Eigen::MatrixXd M(S.rows(), d);
  M.noalias() = S * A;

  const double denom0 = A.col(0).dot(M.col(0));
  if (!std::isfinite(denom0) || denom0 == 0.0) {
    Rf_error("makeVexpC: non-finite or zero denominator for first component");
  }

  vexp(0) = M.col(0).squaredNorm() / denom0;
  cvexp(0) = vexp(0);

  for (int i = 1; i < d; ++i) {
    const int k = i + 1;

    const Eigen::MatrixXd Mi = M.leftCols(k);
    const Eigen::MatrixXd Ai = A.leftCols(k);

    // G = A' (S A) = A' M
    Eigen::MatrixXd G(k, k);
    G.noalias() = Ai.transpose() * Mi;

    // Solve instead of explicit inverse: G^{-1}
    const Eigen::MatrixXd Ginv = G.ldlt().solve(Eigen::MatrixXd::Identity(k, k));

    // Use cyclic trace: tr(M G^{-1} M') = tr(G^{-1} (M' M))
    Eigen::MatrixXd H(k, k);
    H.noalias() = Mi.transpose() * Mi;

    cvexp(i) = (Ginv * H).trace();
    vexp(i) = cvexp(i) - cvexp(i - 1);
  }

  return List::create(Named("vexp") = vexp, Named("cvexp") = cvexp);
}

// [[Rcpp::export]]
double mkvexpfirstCompSC(const Eigen::Map<Eigen::VectorXd>& a,
                         const Eigen::Map<Eigen::MatrixXd>& S) {
  // Returns explained variance of the first component a w.r.t. S.

  const Eigen::VectorXd B = S * a;
  const double denom = a.dot(B);
  if (!std::isfinite(denom) || denom == 0.0) {
    Rf_error("mkvexpfirstCompSC: non-finite or zero denominator");
  }
  return B.squaredNorm() / denom;
}

// [[Rcpp::export]]
List makeVexpSC(const Eigen::Map<Eigen::MatrixXd>& A,
               const Eigen::Map<Eigen::MatrixXd>& S) {
  // Same as makeVexpC, but kept as a separate entry point for API compatibility.

  const int d = A.cols();
  Eigen::VectorXd vexp(d);
  Eigen::VectorXd cvexp(d);

  Eigen::MatrixXd M(S.rows(), d);
  M.noalias() = S * A;

  const double denom0 = A.col(0).dot(M.col(0));
  if (!std::isfinite(denom0) || denom0 == 0.0) {
    Rf_error("makeVexpSC: non-finite or zero denominator for first component");
  }

  vexp(0) = M.col(0).squaredNorm() / denom0;
  cvexp(0) = vexp(0);

  for (int i = 1; i < d; ++i) {
    const int k = i + 1;

    const Eigen::MatrixXd Mi = M.leftCols(k);
    const Eigen::MatrixXd Ai = A.leftCols(k);

    Eigen::MatrixXd G(k, k);
    G.noalias() = Ai.transpose() * Mi;

    const Eigen::MatrixXd Ginv = G.ldlt().solve(Eigen::MatrixXd::Identity(k, k));

    Eigen::MatrixXd H(k, k);
    H.noalias() = Mi.transpose() * Mi;

    cvexp(i) = (Ginv * H).trace();
    vexp(i) = cvexp(i) - cvexp(i - 1);
  }

  return List::create(Named("vexp") = vexp, Named("cvexp") = cvexp);
}

// Not exported (kept consistent with your original file)
List makeEvexpOneCompC(const Eigen::Map<Eigen::MatrixXd>& A,
                       const Eigen::Map<Eigen::MatrixXd>& T,
                       const Eigen::Map<Eigen::MatrixXd>& S,
                       double oldcvexp = 0) {
  // Computes the extra variance explained by one component.

  double cvexp, vexp;

  const int k = T.cols();
  Eigen::MatrixXd G(k, k);
  G.noalias() = T.transpose() * T;

  const Eigen::MatrixXd Ginv = G.ldlt().solve(Eigen::MatrixXd::Identity(k, k));

  // M = T (T'T)^{-1} A' S
  Eigen::MatrixXd AS(A.cols(), S.cols());
  AS.noalias() = A.transpose() * S;

  Eigen::MatrixXd M(T.rows(), AS.cols());
  M.noalias() = T * Ginv * AS;

  cvexp = M.squaredNorm();
  vexp = cvexp - oldcvexp;

  return List::create(Named("vexp") = vexp, Named("cvexp") = cvexp);
}

// [[Rcpp::export]]
Eigen::MatrixXd MkCorCompMat(const Eigen::Map<Eigen::MatrixXd>& A,
                             const Eigen::Map<Eigen::MatrixXd>& S,
                             int d = 0) {
  // Computes the correlation between components (columns of A) w.r.t. S.
  // Returns a d x d matrix.

  if ((d > A.cols()) || (d == 0)) {
    d = A.cols();
  }
  if (d < 2) {
    Rf_error("MkCorCompMat: d must be larger than 1");
  }

  const Eigen::MatrixXd Ad = A.leftCols(d);

  Eigen::MatrixXd SA(S.rows(), d);
  SA.noalias() = S * Ad;

  Eigen::MatrixXd C(d, d);
  C.noalias() = Ad.transpose() * SA;

  const Eigen::VectorXd s = C.diagonal().array().sqrt().matrix();

  C.array().rowwise() /= s.array().transpose();
  C.array().colwise() /= s.array();

  return C;
}

//====================================================
// VIF helpers
//====================================================

static void vifintS(const Eigen::Ref<const Eigen::MatrixXd>& S,
                    int ind, double& vifout) {
  // Computes VIF for variable with index ind (base 1)

  const int ind0 = ind - 1;
  const int p = S.cols();

  if (ind0 > (p - 1)) {
    Rf_error("vifintS: must pass indices less than or equal to p");
  }

  Eigen::MatrixXd T(p - 1, p - 1);
  Eigen::VectorXd s(p - 1), sc(p);

  if (ind0 == 0) {
    T = S.bottomRightCorner(p - 1, p - 1);
    sc = S.col(0);
    s = sc.tail(p - 1);
  }
  if (ind0 == (p - 1)) {
    T = S.topLeftCorner(p - 1, p - 1);
    sc = S.col(p - 1);
    s = sc.head(p - 1);
  }
  if ((ind0 > 0) && (ind0 < (p - 1))) {
    const Eigen::MatrixXd A = S.topLeftCorner(ind0, ind0);
    const Eigen::MatrixXd B = S.topRightCorner(ind0, p - ind0 - 1);
    const Eigen::MatrixXd C = S.bottomRightCorner(p - ind0 - 1, p - ind0 - 1);

    T.topLeftCorner(ind0, ind0) = A;
    T.topRightCorner(ind0, p - ind0 - 1) = B;
    T.bottomLeftCorner(p - ind0 - 1, ind0) = B.transpose();
    T.bottomRightCorner(p - ind0 - 1, p - ind0 - 1) = C;

    sc = S.col(ind0);
    s.head(ind0) = sc.head(ind0);
    s.tail(p - ind0 - 1) = sc.tail(p - ind0 - 1);
  }

  const Eigen::MatrixXd Tm = T.llt().solve(Eigen::MatrixXd::Identity(p - 1, p - 1));
  vifout = (s.transpose() * Tm * s)(0, 0) / sc(ind0);
}

static void vifintSPseudo(const Eigen::Ref<const Eigen::MatrixXd>& S,
                          int ind, double& vifout) {
  // Computes VIF for variable with index ind (base 1) using pseudo-inverse.

  const int ind0 = ind - 1;
  const int p = S.cols();

  if (ind0 > (p - 1)) {
    Rf_error("vifintSPseudo: must pass indices less than or equal to p");
  }

  Eigen::MatrixXd T(p - 1, p - 1);
  Eigen::VectorXd s(p - 1), sc(p);

  if (ind0 == 0) {
    T = S.bottomRightCorner(p - 1, p - 1);
    sc = S.col(0);
    s = sc.tail(p - 1);
  }
  if (ind0 == (p - 1)) {
    T = S.topLeftCorner(p - 1, p - 1);
    sc = S.col(p - 1);
    s = sc.head(p - 1);
  }
  if ((ind0 > 0) && (ind0 < (p - 1))) {
    const Eigen::MatrixXd A = S.topLeftCorner(ind0, ind0);
    const Eigen::MatrixXd B = S.topRightCorner(ind0, p - ind0 - 1);
    const Eigen::MatrixXd C = S.bottomRightCorner(p - ind0 - 1, p - ind0 - 1);

    T.topLeftCorner(ind0, ind0) = A;
    T.topRightCorner(ind0, p - ind0 - 1) = B;
    T.bottomLeftCorner(p - ind0 - 1, ind0) = B.transpose();
    T.bottomRightCorner(p - ind0 - 1, p - ind0 - 1) = C;

    sc = S.col(ind0);
    s.head(ind0) = sc.head(ind0);
    s.tail(p - ind0 - 1) = sc.tail(p - ind0 - 1);
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);

  int rank = 0;
  const Eigen::VectorXd la = es.eigenvalues();
  Eigen::VectorXd lai(p - 1);
  lai.setZero();

  for (int i = 0; i < p - 1; ++i) {
    if (la(i) > 1e-5) {
      rank += 1;
      lai(p - i - 2) = 1.0 / la(i);
    }
  }

  const Eigen::MatrixXd V = es.eigenvectors().rightCols(rank);
  const Eigen::MatrixXd Tm = V * lai.head(rank).reverse().asDiagonal() * V.transpose();

  vifout = (s.transpose() * Tm * s)(0, 0) / sc(ind0);
}

// [[Rcpp::export]]
Eigen::VectorXd vifSC(const Eigen::Map<Eigen::MatrixXd>& S,
                      const Eigen::VectorXi& ind) {
  // Computes VIF for all variables listed in ind using covariance/correlation S.
  // ind is expected base 1.

  const int p = S.cols();
  const int n = S.rows();

  if (p == 1) {
    Rf_error("vifSC: must pass a matrix with at least 2 columns");
  }
  if (p != n) {
    Rf_error("vifSC: must pass a variance or correlation matrix");
  }

  const int nind = ind.size();

  if ((nind > p) || (ind.maxCoeff() > p)) {
    Rf_error("vifSC: indices must be base 1 and <= p");
  }

  Eigen::VectorXd vif(nind);
  for (int i = 0; i < nind; ++i) {
    vifintS(S, ind(i), vif(i));
  }

  return vif;
}

// [[Rcpp::export]]
Eigen::VectorXd vifSPseudoC(const Eigen::Map<Eigen::MatrixXd>& S,
                            const Eigen::VectorXi& ind) {
  // Computes VIF for all variables listed in ind using covariance/correlation S.
  // ind is expected base 1.

  const int p = S.cols();
  const int n = S.rows();

  if (p == 1) {
    Rf_error("vifSPseudoC: must pass a matrix with at least 2 columns");
  }
  if (p != n) {
    Rf_error("vifSPseudoC: must pass a variance or correlation matrix");
  }

  const int nind = ind.size();

  if ((nind > p) || (ind.maxCoeff() > p)) {
    Rf_error("vifSPseudoC: indices must be base 1 and <= p");
  }

  Eigen::VectorXd vif(nind);
  for (int i = 0; i < nind; ++i) {
    vifintSPseudo(S, ind(i), vif(i));
  }

  return vif;
}
