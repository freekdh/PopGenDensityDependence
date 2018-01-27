#include <Rcpp.h>
#include <vector>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/binomial.hpp>

enum parnames {n, freq, prob};

struct pars{
  double z1,z2,k;
};

std::vector<std::array<double, 3>>::iterator it;

double s(const double &n, const double &p, const pars &par){
  return (((par.k)- n * p - n * (1 - p)) * (par.z1 - par.z2))/((par.k) + ((par.k) - n*p-n*(1-p))*(p*(par.z1 - par.z2)+ par.z2));
}


double n2(const double &n, const double &p, const pars &par){
  if(n==0.0) return 0.0;
  return static_cast<unsigned int>(n + n*((par.z1)*p+(par.z2)*(1.0-p))*(1.0-(n/(par.k))));
}

  
void ItteratorNP(std::vector<std::array<double, 3> > &input, std::vector<std::array<double, 3> > &output, const pars &par){
  output.clear();
  for(it = input.begin(); it != input.end(); ++it){
    std::array<double,3> pre, nex;
    pre[n] = (*it)[n];
    pre[freq] = (*it)[freq];
    pre[prob] = (*it)[prob];
    
    nex[n] = n2(pre[n],pre[freq],par);
    for(unsigned int i = 0; i < nex[n]; ++i){
      double snow = s(pre[n], pre[freq], par);
      nex[freq] = (double)i / nex[0];
      nex[prob] = boost::math::binomial_coefficient<double>((int)nex[n], (double)i) * 
        pow((-1.0 + pre[freq]) * (-1.0 + pre[freq] * snow), -i + nex[n]) * 
        pow((pre[freq] - (-1.0 + pre[freq]) * pre[freq] * snow), (double)i);

      output.push_back(nex);
    }
  }
}

// [[Rcpp::export]]
NumericVector Solve(double n1, double n2, double z1, double z2, double k) {
  const struct pars par = {z1, z2, k};
  
  std::vector<std::array<double,3 > > input, output;
  std::array<double, 3> trial = {0.1,0.2, 0.4};
  
  input.push_back(trial);  

  ItteratorNP(input, output, par);
  
  return 0;
}

