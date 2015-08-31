#ifndef __PGS_TOOLS_INCLUDED__
#define __PGS_TOOLS_INCLUDED__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

vec countRep_cpp(vec);

List indGen_cpp(vec);

uvec seqJoin_vec(uvec, uvec, vec);
  
uvec seqJoin_int(uvec, int);
  
mat corr_est_normal_cpp(vec, mat, vec, vec, int, std::string);
    
#endif