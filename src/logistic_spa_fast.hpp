#ifndef INCLUDE_LOGISTIC
#define INCLUDE_LOGISTIC

#include <RcppArmadillo.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <sstream>
#include <time.h>
#include <Rcpp.h>
#include <stdint.h>
#include <boost/math/distributions/normal.hpp>

double Korg_Binom(double t1, arma::vec &mu, arma::vec &g);
double K1_adj_Binom(double t1, arma::vec &mu, arma::vec &g, double q);
double K2_Binom(double t1, arma::vec &mu, arma::vec &g);
Rcpp::List getroot_K1_Binom(double init, arma::vec &mu, arma::vec &g, double q, double tol, int maxiter = 1000);
Rcpp::List Get_Saddle_Prob_Binom(double zeta, arma::vec &mu, arma::vec &g, double q, bool logp = false);
Rcpp::List SPA_binary(arma::vec &mu, arma::vec &g, double q, double qinv, double pval_noadj, double tol, bool logp = false);
double sum_arma1(arma::vec &X);
double add_logp(double p1, double p2);

#endif
