#include <Rcpp.h>
#include <vector>
#include <assert.h>

using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/binomial.hpp>

// [[Rcpp::export]]
inline double s(const double &n, const double &p, const double &k, const double &z1, const double &z2){
  return ((k - n)*(z1 - z2))/(k + k*p*z1 - n*p*z1 - (k-n)*(-1.0+p)*z2);
}

// [[Rcpp::export]]
inline int n2(const double &n, const double &p, const double &k, const double &z1, const double &z2){
  if(n==0.0) return 0.0;
  
  return static_cast<int>(n + n*((z1)*p+(z2)*(1.0-p))*(1.0-(n/(k))));

}

// [[Rcpp::export]]
NumericMatrix ItteratorNP(NumericMatrix input, const double &z1, const double &z2, const int &ngen ){
  
  int k = input.nrow();
  NumericMatrix output(k,k);
  
  std::vector<std::vector<double>> invec(k, std::vector<double>(k));
  std::vector<std::vector<double>> outvec(k, std::vector<double>(k,0.0));
  
  
  for(int i = 0; i < k; ++i){
    for(int j = 0; j < k; ++j){
      invec[i][j] = input(i,j);
    }
  }
  
  for(int ngenit = 0 ; ngenit < ngen; ++ngenit){
    for(int i = 0; i < k; ++i){
      for(int j = 0; j < k; ++j){
        if(invec[i][j] == 0.0) {continue;}
        else{
          int nnext = n2(i+j,(double)i/(double)(i+j),k, z1, z2);
          double f = (double)i/(double)(i+j);
          double snow = s(i+j, f, k, z1, z2);
          for(int o = 0; o <= nnext; ++o){
            outvec[o][nnext-o] += invec[i][j] * ((boost::math::binomial_coefficient<double>(nnext, o)) * 
              pow(((-1.0 + f)*(-1.0 + f*snow)),(-o + nnext))* pow((f - (-1.0 + f)*f*snow),o));
          }
        }
      }
    }
    
    for(int i = 0; i < k; ++i){
      for(int j = 0; j < k; ++j){
        invec[i][j] = outvec[i][j];
        outvec[i][j] = 0.0;
      }
    }   
  }
  
  
  for(int i = 0; i < k; ++i){
    for(int j = 0; j < k; ++j){
      output(i,j) = invec[i][j];
    }
  }
  
    
  return output;
}

/*
// [[Rcpp::export]]
NumericMatrix Solve(const int n1, const int n2, const double z1, const double z2, const double k) {

  std::vector<std::vector<double> > input(k, std::vector<double>(k,0.0));
  std::vector<std::vector<double> > output(k, std::vector<double>(k,0.0));
  
  input[n1][n2] = 1.0;
    
  ItteratorNP(input, output, k, z1, z2);
  
  NumericMatrix outR(k, k);
  
  for(int i = 0; i < k; ++i){
    for(int j = 0; j < k; ++j){
      outR(i,j) = output[i][j];
    }
  }
  
  return outR;
}
 */

