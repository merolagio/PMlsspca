#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <iostream>


using namespace Rcpp;
using namespace Eigen;
using namespace std;

//dont forget
//Rcpp::compileAttributes()


// [[Rcpp::depends(RcppEigen)]]

//
using Eigen::Map;               	// 'maps' rather than copies
using Eigen::Matrix;                  //  matrix generic
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::Transpositions;
using Eigen::HouseholderQR;    // Fast scalable QR solver
using Eigen::ColPivHouseholderQR;    // Fast scalable QR solver
using Eigen::FullPivHouseholderQR; // slow full (colsand rows pivoting)
using Eigen::GeneralizedSelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::JacobiSVD;
using Eigen::LLT;
using Eigen::LDLT;
using Rcpp::List;
using Rcpp::wrap;


// new ===========================================

// // //' @export
 // [[Rcpp::export]]  
 List deflXCforR(Eigen::VectorXd comp, Eigen::MatrixXd K){
   // # deflates data matrix K of comp
   // ## ===================================================
   // comp component
   // K = matrix (possibly deflated of previous comps)
   // #  
   // K <-- {}(K - t*t'/(t't)}K // deflated K matrix
   // ## ===================================================
   
   double vexp = 0;
   if( K.rows() != comp.size())
     Rf_error("matrix and component are not comptible for deflation");
   
   double tt = comp.squaredNorm();
   Eigen::MatrixXd T = (comp*comp.transpose()).array()/tt;
   Eigen::MatrixXd  M = T*K;
   vexp = M.squaredNorm();
   K = K - M;
   
   return List::create (Named("K") = K, Named("vexp") = vexp);
 }

// [[Rcpp::export]]
List makeVexpC(Eigen::MatrixXd A, Eigen::MatrixXd S){
  // computes the extra variance explained by a set of components
  int p = A.cols();
  Eigen::VectorXd vexp(p);
  Eigen::VectorXd cvexp(p);
  Eigen::MatrixXd M = S * A; 
  
  // Rcout << "sqrd norm Sa_1\n" <<  (M.col(0)).squaredNorm() << endl;
  // Rcout << "sqrd norm a_1\n" <<  A.col(0).squaredNorm() << endl;
  vexp(0) = M.col(0).squaredNorm()/(A.col(0).transpose() * M.col(0));
  cvexp(0) = vexp(0);
  for (int i = 1; i < p; i++){
    cvexp(i) = (M.leftCols(i + 1) *  (A.leftCols(i + 1).transpose() * M.leftCols(i + 1)).inverse() *
      M.leftCols(i + 1).transpose()).trace();
    vexp(i) = cvexp(i) - cvexp(i -1);
  }
  return List::create(Named("vexp") = vexp, Named("cvexp") = cvexp );
}

// [[Rcpp::export]]
 double mkvexpfirstCompSC(Eigen::VectorXd a, Eigen::MatrixXd S){
   // returns the variance explained by the first compenent
   Eigen::MatrixXd B = S * a;
   
   return((B.transpose() * B)(0,0)/(a.transpose() * B)(0,0));
 }


// // [[Rcpp::export]]
List makeEvexpOneCompC(Eigen::MatrixXd A, Eigen::MatrixXd T, Eigen::MatrixXd S, 
                       double oldcvexp = 0){
  
  // computes the extra variance explained by one components
  
  double cvexp, vexp;
  Eigen::MatrixXd M = T * (T.transpose() * T).inverse() * A.transpose() * S; 
  cvexp = M.squaredNorm();
  // Rcout << "sqrd norm Sa_1\n" <<  (M.col(0)).squaredNorm() << endl;
  // Rcout << "sqrd norm a_1\n" <<  A.col(0).squaredNorm() << endl;
  vexp = cvexp - oldcvexp;
  return List::create(Named("vexp") = vexp, Named("cvexp") = cvexp );
}

//computes the correlation between components
//returns Cor Matrix
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd MkCorCompMat(Eigen::MatrixXd A, Eigen::MatrixXd S, int d = 0){
   if ((d > A.cols()) | (d == 0))
     d = A.cols(); 
   if (d < 2)
     Rf_error("MkCorCompMat: d must be larger than 1");
   
   Eigen::MatrixXd C(d, d);
   C = A.leftCols(d).transpose() * S * A.leftCols(d);
   Eigen::VectorXd s(C.diagonal().array().sqrt());
   
   C.array().rowwise() /= s.array().transpose();
   C.array().colwise() /= s.array();
   return C;
 }
//end new

//=========================================
//               Basic matrix operations
//=============================================

// a, b, c for matrices
// u and v for vectors
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd aatC(Eigen::Map<Eigen::MatrixXd > A){
   int p(A.rows());
   Eigen::MatrixXd M(p,p);
   M.noalias() = A * A.adjoint();  
   return M;  
 }
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd ataC(Eigen::Map<Eigen::MatrixXd > A){
   int p(A.rows());
   Eigen::MatrixXd M(p,p);
   M.noalias() = A.adjoint() * A;  
   return M;  
 }
// At S
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd atbC(Eigen::Map<Eigen::MatrixXd > A, Eigen::Map<Eigen::MatrixXd > B){
   int p(A.cols()), q(B.cols());
   Eigen::MatrixXd M(p, q);
   M.noalias() =  A.adjoint() * B;  
   return M;  
 }
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd abC(Eigen::Map<Eigen::MatrixXd > A, Eigen::Map<Eigen::MatrixXd > B){
   int p(A.rows());
   int q(B.cols());
   int n1(A.cols());
   int n2(B.rows());
   if ( n1 != n2)
     Rf_error("abc: A and B must be compatible for multiplication A * B");
   Eigen::MatrixXd M(p,q);
   M.noalias() = A * B;  
   return M;  
 }

//A t(B)
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd abtC(Eigen::Map<Eigen::MatrixXd > A, Eigen::Map<Eigen::MatrixXd > B){
   int p(A.rows()), q(B.rows());
   if(A.cols() != B.cols())
     Rf_error(" abtC: A and B must have the same number of columns");
   Eigen::MatrixXd M(p,q);
   M.noalias() = A * B.adjoint();  
   return M;  
 }

//t(A) D A
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd atdaC(Eigen::Map<Eigen::MatrixXd > A, Eigen::Map<Eigen::MatrixXd > D){
   int p(A.cols());
   Eigen::MatrixXd M(p,p);
   M.noalias() =  A.adjoint() * D * A;  
   return M;  
 }

//t(A) D B
// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd atdbC(Eigen::Map<Eigen::MatrixXd > A, Eigen::Map<Eigen::MatrixXd > D,
                       Eigen::Map<Eigen::MatrixXd > B){
   int p(A.cols());
   Eigen::MatrixXd M(p,p);
   M.noalias() =  A.adjoint() * D * B;  
   return M;  
 }


// nome sbagliato, dovrebbe essere AvC
// //' @export
 // [[Rcpp::export]]  
 Eigen::VectorXd avC(Eigen::MatrixXd A, Eigen::VectorXd v){
   
   return((A * v).col(0));
 }

// nome sbagliato, dovrebbe essere vtAC
// //' @export
 // [[Rcpp::export]]  
 Eigen::VectorXd vtaC(Eigen::MatrixXd A, Eigen::VectorXd v){
   
   return((v.transpose() * A ).row(0));
 }


// //' @export
 // [[Rcpp::export]]  
 double vtvC(Eigen::VectorXd v){
   
   return(v.squaredNorm());
 }

// //' @export
 // [[Rcpp::export]]  
 double vtuC(Eigen::VectorXd v, Eigen::VectorXd u){
   
   return((v.transpose() * u)(0,0));
 }

// //' @export
 // [[Rcpp::export]]
 Eigen::ArrayXXd scaleC(Eigen::Map<Eigen::MatrixXd> A, bool center = true, bool scale = true){
   
   Eigen::MatrixXd C(A);
   Eigen::VectorXd m(A.cols());
   if (center){
     m = C.colwise().mean();
     C.array().rowwise() -= m.array().transpose();//
   }  
   if (scale){
     m = C.colwise().norm();
     C.array().rowwise() /= m.array().transpose();
   }
   return C;
 }

// //' @export
 // [[Rcpp::export]]
 Eigen::MatrixXd scaleColsC(Eigen::Map<Eigen::MatrixXd > A, int normtype,
                            Eigen::Map<Eigen::VectorXd > sig){ 
   //  sig is vector of coeff to multiply each column (usually +/- 1)
   // non serve .array()
   Eigen::ArrayXXd M(A);
   if (normtype == 2){
     Eigen::VectorXd l2 = A.colwise().norm();//.array()
     //    std::cout << l2 << std::endl;
     M.rowwise() /= (l2.array() * sig.array()).transpose() ;
   }
   else{
     if (normtype == 1){
       Eigen::VectorXd l1 = A.array().abs().colwise().sum();
       M.rowwise() /= (l1.array() * sig.array()).transpose();//
     }
     else {
       M.rowwise() /= sig.array().transpose();
     }
   }  
   return M;
 }  

// //' @export
 // [[Rcpp::export]]
 NumericVector sortC(NumericVector v) {
   NumericVector y = clone(v);
   std::sort(y.begin(), y.end());
   return y;
 }

// //' @export
 // [[Rcpp::export]]
 double sumdiagC(Eigen::Map<Eigen::MatrixXd > S){
   return  S.trace();
 }


//======================================
//        helpers for lsspca
//======================================

// //' @export 
 // [[Rcpp::export]]
 List EigenC(Eigen::Map<Eigen::MatrixXd> S) {
   if (S.cols() != S.rows())
     Rcpp::stop("S in EigenC must be symmetric");
   SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
   return List::create(Named("vec") = Eigen::MatrixXd(es.eigenvectors()).rowwise().reverse(),
                       Named("val") = Eigen::VectorXd(es.eigenvalues()).reverse());
 }

// EigenC computes generalised eigendecomp A x = Bx mu
// //' @export 
 // [[Rcpp::export]]
 List GenEigenC(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B) {
   GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es(A,B);
   return List::create(Named("vec") = Eigen::MatrixXd(es.eigenvectors()).rowwise().reverse(),
                       Named("val") = Eigen::VectorXd(es.eigenvalues()).reverse());
 }

// SolveC computes inverse
// //' @export 
 // [[Rcpp::export]]
 Eigen::MatrixXd solveC(Eigen::MatrixXd S) {
   const int d = S.cols();
   const Eigen::MatrixXd A(S.llt().solve(MatrixXd::Identity(d,d)));
   return A;
 }

// utilities========================================


void vifintS(Eigen::MatrixXd S, int ind, double& vifout){
  // computes vif for var with index ind 
  // takes indices base 1, converts them to base 0 
  int ind0 = ind - 1;
  int p = S.cols();
  if (ind0 > (p - 1))
    Rf_error("VifintS: must pass indices less than p");
  Eigen::MatrixXd T(p - 1, p - 1);
  Eigen::VectorXd s(p - 1), sc(p);
  if (ind0 == 0){
    T = S.bottomRightCorner(p - 1, p - 1);
    sc = S.col(0);
    s = sc.tail(p - 1);
  }
  if (ind0 == (p - 1)){
    T = S.topLeftCorner(p - 1, p - 1);
    sc = S.col(p - 1);
    s = sc.head(p - 1);
  }
  if((ind0 > 0) & (ind0 < (p - 1))){
    Eigen::MatrixXd A = S.topLeftCorner(ind0, ind0);
    Eigen::MatrixXd B = S.topRightCorner(ind0, p - ind0 - 1);
    Eigen::MatrixXd C = S.bottomRightCorner(p - ind0 - 1, p - ind0 - 1);
    
    T.topLeftCorner(ind0, ind0) = A;
    T.topRightCorner(ind0, p - ind0 - 1) = B;
    T.bottomLeftCorner(p - ind0 - 1, ind0) = B.transpose();
    T.bottomRightCorner(p - ind0 - 1, p - ind0 - 1) = C;
    sc = S.col(ind0);
    s.head(ind0) = sc.head(ind0);
    s.tail(p - ind0 - 1) = sc.tail(p - ind0 - 1);
  }
  const Eigen::MatrixXd Tm(T.llt().solve(MatrixXd::Identity(p - 1, p - 1)));
  
  vifout = (s.transpose() * Tm * s)(0, 0)/sc(ind0);
}

//  // ' @export
//  // [[Rcpp::export]]
void vifintSPseudo(Eigen::MatrixXd S, int ind, double& vifout){
  // computes vif for var with index ind 
  // takes indices base 1, converts them to base 0 
  
  int ind0 = ind - 1;
  int p = S.cols();
  if (ind0 > (p - 1))
    Rf_error("VifintS: must pass indices less than p");
  Eigen::MatrixXd T(p - 1, p - 1);
  Eigen::VectorXd s(p - 1), sc(p);
  if (ind0 == 0){
    T = S.bottomRightCorner(p - 1, p - 1);
    sc = S.col(0);
    s = sc.tail(p - 1);
  }
  if (ind0 == (p - 1)){
    T = S.topLeftCorner(p - 1, p - 1);
    sc = S.col(p - 1);
    s = sc.head(p - 1);
  }
  if((ind0 > 0) & (ind0 < (p - 1))){
    Eigen::MatrixXd A = S.topLeftCorner(ind0, ind0);
    Eigen::MatrixXd B = S.topRightCorner(ind0, p - ind0 - 1);
    Eigen::MatrixXd C = S.bottomRightCorner(p - ind0 - 1, p - ind0 - 1);
    
    T.topLeftCorner(ind0, ind0) = A;
    T.topRightCorner(ind0, p - ind0 - 1) = B;
    T.bottomLeftCorner(p - ind0 - 1, ind0) = B.transpose();
    T.bottomRightCorner(p - ind0 - 1, p - ind0 - 1) = C;
    sc = S.col(ind0);
    s.head(ind0) = sc.head(ind0);
    s.tail(p - ind0 - 1) = sc.tail(p - ind0 - 1);
  }
  //  const Eigen::MatrixXd Tm(T.llt().solve(MatrixXd::Identity(p - 1, p - 1)));
  
  // does pseudo inverse with eigen :-(
  // not sure does pseudo inv   const Eigen::MatrixXd Tm(T.inverse());
  SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
  
  // remember size S is p - 1
  int rank = 0;
  Eigen::VectorXd la = es.eigenvalues();
  Eigen::VectorXd lai(p - 1);
  for(int i = 0; i < p - 1; ++i){
    if(la(i) > 1e-5){
      rank = rank + 1;
      lai(p - i - 2) = 1/la(i);
      // Rcout << " p - i - 2 = " << p - i - 2 << endl;
    }
  }
  // Rcout << "rank is " << rank << ";  p - rank = " << p - rank << endl;
  
  Eigen::MatrixXd V = es.eigenvectors().rightCols(rank);
  
  // Rcout << " eigenvalues " << la.tail(rank).transpose() << endl << endl;
  // Rcout << " inv eigenvalues " << lai.head(rank).reverse().transpose() << endl;
  
  //   Eigen::VectorXd lais = es.eigenvalues();
  const Eigen::MatrixXd Tm = V * lai.head(rank).reverse().asDiagonal() * V.transpose(); 
  vifout = (s.transpose() * Tm * s)(0, 0)/sc(ind0);
}

// //' @export
 // [[Rcpp::export]]
 Eigen::VectorXd vifSC(Eigen::Map<Eigen::MatrixXd> S, Eigen::VectorXi ind){
   // computes vif for all variables in S, cov or cor matrix
   // takes indices base 1, converted to base 0 in internal
   int p = S.cols(), n = S.rows();
   
   if (p == 1){
     Rf_error("must pass a matrix with at least 2 columns to vifSC");
   }
   if (p != n){
     Rf_error("must pass a variance or correlation matrix to vifSC");
   }
   
   int nind = ind.size();
   
   if ((nind > p ) | (ind.maxCoeff() > p)){
     Rf_error("something wrong with indices for vifC, it takes indices base 0");
   }
   
   Eigen::VectorXd vif(nind);
   for (int i = 0; i < nind; ++i){
     vifintS(S, ind(i), vif(i));
   }
   return(vif);
 }
// //' @export
 // [[Rcpp::export]]
 Eigen::VectorXd vifSPseudoC(Eigen::Map<Eigen::MatrixXd> S, Eigen::VectorXi ind){
   // computes vif for all variables in S, cov or cor matrix
   // takes indices base 1, converted to base 0 in internal
   int p = S.cols(), n = S.rows();
   
   if (p == 1){
     Rf_error("must pass a matrix with at least 2 columns to vifSC");
   }
   if (p != n){
     Rf_error("must pass a variance or correlation matrix to vifSC");
   }
   
   int nind = ind.size();
   
   if ((nind > p ) | (ind.maxCoeff() > p)){
     Rf_error("something wrong with indices for vifC, it takes indices base 0");
   }
   
   Eigen::VectorXd vif(nind);
   for (int i = 0; i < nind; ++i){
     vifintSPseudo(S, ind(i), vif(i));
   }
   return(vif);
 }


