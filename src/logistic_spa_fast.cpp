// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include "logistic_spa_fast.hpp"
// boost::math::normal norm;

// Compute p value

// [[Rcpp::export]]
Rcpp::List SPA_binary(arma::vec &mu, arma::vec &g, double q, double qinv, double pval_noadj, double tol, int maxiter, bool logp)
{
  using namespace Rcpp;
  List result;
  double p1, p2, pval;
  bool Isconverge = true;
  Rcpp::List outuni1 = getroot_K1_Binom(0, mu, g, q, tol, maxiter);
  Rcpp::List outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol, maxiter);
  // double outuni1root = outuni1["root"];
  // double outuni2root = outuni2["root"];
  // bool Isconverge1 = outuni1["Isconverge"];
  // bool Isconverge2 = outuni2["Isconverge"];

  // std::cout << "outuni1root" << outuni1root << std::endl;
  // std::cout << "outuni2root" << outuni2root << std::endl;
  // std::cout << "Isconverge1" << Isconverge1 << std::endl;
  // std::cout << "Isconverge2" << Isconverge2 << std::endl;

  Rcpp::List getSaddle;
  Rcpp::List getSaddle2;
  if (outuni1["Isconverge"] && outuni2["Isconverge"])
  {
    // std::cout << "q is " << q << " 3 qinv is " << qinv << std::endl;
    getSaddle = Get_Saddle_Prob_Binom(outuni1["root"], mu, g, q, logp);

    if (getSaddle["isSaddle"])
    {
      p1 = getSaddle["pval"];
    }
    else
    {

      if (logp)
      {
        p1 = pval_noadj - std::log(2);
      }
      else
      {
        p1 = pval_noadj / 2;
      }
    }
    getSaddle2 = Get_Saddle_Prob_Binom(outuni2["root"], mu, g, qinv, logp);
    if (getSaddle2["isSaddle"])
    {
      p2 = getSaddle2["pval"];
    }
    else
    {
      if (logp)
      {
        p2 = pval_noadj - std::log(2);
      }
      else
      {
        p2 = pval_noadj / 2;
      }
    }

    if (logp)
    {
      pval = add_logp(p1, p2);
    }
    else
    {
      pval = std::abs(p1) + std::abs(p2);
    }
    Isconverge = true;
  }
  else
  {
    pval = pval_noadj;
    Isconverge = false;
  }
  result["pvalue"] = pval;
  result["Isconverge"] = Isconverge;
  return (result);
}

// Compute CGF function

double Korg_Binom(double t1, arma::vec &mu, arma::vec &g)
{
  arma::vec temp = arma::log(1 - mu + mu % arma::exp(g * t1));
  double out = arma::sum(temp);

  return (out);
}

double K1_adj_Binom(double t1, arma::vec &mu, arma::vec &g, double q)
{
  arma::vec temp1;
  arma::vec temp2;
  arma::vec temp3;

  temp1 = (1 - mu) % arma::exp(-g * t1) + mu;
  temp2 = mu % g;
  temp3 = temp2 / temp1;
  double out = arma::sum(temp3) - q;

  return (out);
}

double K2_Binom(double t1, arma::vec &mu, arma::vec &g)
{
  arma::vec temp0;
  arma::vec temp1;
  arma::vec temp2;
  arma::vec temp3;

  temp0 = arma::exp(-g * t1);
  temp1 = (1 - mu) % temp0 + mu;
  temp1 = arma::pow(temp1, 2);
  temp2 = arma::pow(g, 2) % temp0;
  temp2 = (1 - mu) % mu % temp2;
  temp3 = temp2 / temp1;
  double out = sum_arma1(temp3);

  return (out);
}

Rcpp::List getroot_K1_Binom(double init, arma::vec &mu, arma::vec &g, double q, double tol, int maxiter)
{
  Rcpp::List result;
  double root;
  int niter;
  bool Isconverge;
  int rep;
  double K1_eval, K2_eval, t, tnew, newK1;
  double prevJump;
  double gpos = arma::accu(g.elem(find(g > 0)));
  double gneg = arma::accu(g.elem(find(g < 0)));
  // std::cout << "q is " << q << std::endl;
  // std::cout << "gpos is " << gpos << std::endl;
  // std::cout << "gneg is " << gneg << std::endl;
  if (q >= gpos || q <= gneg)
  {
    root = std::numeric_limits<double>::infinity();
    niter = 0;
    Isconverge = true;
  }
  else
  {
    t = init;
    K1_eval = K1_adj_Binom(t, mu, g, q);
    prevJump = std::numeric_limits<double>::infinity();
    int rep = 1;
    bool conv = true;
    while (rep <= maxiter)
    {
      K2_eval = K2_Binom(t, mu, g);
      tnew = t - K1_eval / K2_eval;
      if (tnew == NA_REAL)
      {
        conv = false;
        break;
      }

      if (std::abs(tnew - t) < tol)
      {
        conv = true;
        break;
      }

      if (rep == maxiter)
      {
        conv = false;
        break;
      }

      newK1 = K1_adj_Binom(tnew, mu, g, q);
      if (arma::sign(K1_eval) != arma::sign(newK1))
      {
        if (std::abs(tnew - t) > (prevJump - tol))
        {
          tnew = t + (arma::sign(newK1 - K1_eval)) * prevJump / 2;
          newK1 = K1_adj_Binom(tnew, mu, g, q);
          prevJump = prevJump / 2;
        }
        else
        {
          prevJump = std::abs(tnew - t);
        }
      }

      rep = rep + 1;
      t = tnew;
      K1_eval = newK1;
    }
    root = t;
    niter = rep;
    Isconverge = conv;
  }
  // std::cout << "here13" << std::endl;
  // std::cout << "root: " << root << std::endl;
  // std::cout << "niter: " << niter << std::endl;
  result["root"] = root;
  result["niter"] = niter;
  result["Isconverge"] = Isconverge;
  return (result);
}

Rcpp::List Get_Saddle_Prob_Binom(double zeta, arma::vec &mu, arma::vec &g, double q, bool logp)
{
  double k1 = Korg_Binom(zeta, mu, g);
  double k2 = K2_Binom(zeta, mu, g);
  double temp1, w, v, Ztest, pval;
  double negative_infinity = -std::numeric_limits<double>::infinity();

  temp1 = zeta * q - k1;
  Rcpp::List result;
  bool isSaddle = false;
  // std::cout << "k1 " << k1 << std::endl;
  // std::cout << "k2 " << k2 << std::endl;
  // std::cout << "temp1 " << temp1 << std::endl;
  // std::cout << "zeta " << zeta << std::endl;
  // std::cout << "q " << q << std::endl;

  bool flagrun = false;
  if (std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0)
  {
    w = arma::sign(zeta) * std::sqrt(2 * temp1);
    v = zeta * std::sqrt(k2);
    if (w != 0)
    {
      flagrun = true;
    }
  }

  // w = arma::sign(zeta) * std::sqrt(2 *temp1);
  // v = zeta *  std::sqrt(k2);

  // if(std::isfinite(k1) && std::isfinite(k2) && temp1 > 0 && w != 0 && (arma::sign(v) == arma::sign(w)))
  if (flagrun)
  {
    // temp1 = zeta * q - k1;
    Ztest = w + (1 / w) * std::log(v / w);

    // std::cout << "Ztest: " << Ztest << std::endl;
    // std::cout << "w: " << w << std::endl;
    // std::cout << "v: " << v << std::endl;
    // boost::math::normal norm_dist(0,1);

    double pval0;
    if (Ztest > 0)
    {
      // pval0 = boost::math::cdf(complement(norm_dist,Ztest));
      pval0 = R::pnorm(Ztest, 0, 1, false, logp);
      pval = pval0;
      // if(logp){
      //       pval = std::log(pval);
      // }
    }
    else
    {
      // pval0 = boost::math::cdf(norm_dist,Ztest);
      pval0 = R::pnorm(Ztest, 0, 1, true, logp);
      pval = -pval0;

      //      -Rcpp::pnorm( Ztest, mean = 0.0, sd = 1.0, lower = true, log = logp );
    }

    isSaddle = true;
  }
  else
  {
    if (logp)
    {
      pval = negative_infinity;
    }
    else
    {
      pval = 0;
    }
  }
  result["pval"] = pval;
  result["isSaddle"] = isSaddle;
  return (result);
}

// Utilities

double sum_arma1(arma::vec &X)
{
  double sum = 0;
  for (uint i = 0; i < X.n_elem; ++i)
  {
    if (arma::is_finite(X(i)))
      sum += X(i);
  }
  return sum;
}

double add_logp(double p1, double p2)
{
  using namespace Rcpp;
  p1 = -std::abs(p1);
  p2 = -std::abs(p2);
  double maxp = std::max(p1, p2);
  double minp = std::min(p1, p2);
  double result = maxp + std::log(1 + std::exp(minp - maxp));
  return (result);
}
