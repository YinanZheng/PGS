#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "PGS_tools.h"

using namespace arma;
using namespace Rcpp;

//-----------------------------------------------------
// SCAD
vec q_scad_cpp(vec theta, double lambda, double a) {
  
  int p = theta.n_elem;
  vec b1(p), b2(p);
  b1.fill(0); b2.fill(0);
  theta = abs(theta);
  b1.elem(find(theta>lambda)).fill(1);
  b2.elem(find(theta<(lambda * a))).fill(1);
  
  return lambda * (1 - b1) + ((lambda * a) - theta) % b2 % (1 / (a-1) * b1);
}

//-----------------------------------------------------
// S, H, E matrix estimation
List S_H_E_normal_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R_full, vec beta_val, int n, vec m, int obs_n, int p, uvec start, uvec end, vec repm, uvec repm_times, double lam){
  
  double sigma2;
  int i = 0, m_temp, repm_n = repm.n_elem;
  vec Y;
  mat resi, S(p,1), H(p,p), E, X, Xh, hat_R;
  S.fill(0.0); H.fill(0.0);
  
  resi = y_vect - x_mat * beta_val;
  sigma2 = as<double>(wrap(sum(pow(resi, 2))/(obs_n - p)));
  
  for(int j=0; j<repm_n; j++)
  {
    m_temp = repm_times(j)-1;
    hat_R = hat_R_full.submat(0, 0, m_temp, m_temp);  //Trim full hat_R to required m
    for(int k=0; k<repm(j); k++)
    {
      X = x_mat.rows(start(i) , end(i));
      Y = y_vect.subvec(start(i) , end(i));
      Xh = X.t() * hat_R.i();
      S = S + 1 / sigma2 * (Xh * (Y - X * beta_val));
      H = H + 1 / sigma2 * (Xh * X);
      i++;
    }
  }
  
  E = diagmat(q_scad_cpp(beta_val, lam, 3.7) / (abs(beta_val) + 1e-6));
  
  return List::create(Named("S") = S,
                      Named("H") = H,
                      Named("E") = E,
                      Named("hat.R") = hat_R,
                      Named("sigma2") = sigma2);
}

//-----------------------------------------------------
// Shinked beta estimation for a given lambda
List beta_shrink_normal_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R_full, vec beta_ini, int n, vec m, int obs_n, int p, uvec start, uvec end, vec repm, uvec repm_times, double lam, double eps_stop, int iter_try){
  
  int i = 0, step_up = 0, repm_n = repm.n_elem, m_temp;
  vec beta_new = beta_ini, Y, var_sand;
  mat S, H, E, M(p,p), X, EPS, hee, hnee, hat_R;
  M.fill(0.0);
  double flag_stop = math::inf();
  double sigma2;
  
  List S_H_E_val = S_H_E_normal_cpp(y_vect, x_mat, id_vect, hat_R_full, beta_new, n, m, obs_n, p, start, end, repm, repm_times, lam);
  S = as<mat>(S_H_E_val[0]);
  H = as<mat>(S_H_E_val[1]);
  E = as<mat>(S_H_E_val[2]);
  
  while(flag_stop > eps_stop && step_up < iter_try) // When flag_stop < eps_stop, converge, or get out of loop when > max try
  {
    vec beta_old = beta_new;
    beta_new = beta_old + (H + n * E + (1e-6) * eye(p,p)).i() * (S - n * E * beta_old);
    S_H_E_val = S_H_E_normal_cpp(y_vect, x_mat, id_vect, hat_R_full, beta_new, n, m, obs_n, p, start, end, repm, repm_times, lam);
    S = as<mat>(S_H_E_val[0]);
    H = as<mat>(S_H_E_val[1]);
    E = as<mat>(S_H_E_val[2]);
    step_up++;
    flag_stop = sum(abs(beta_old-beta_new));
  }
  hat_R_full = as<mat>(S_H_E_val[3]);
  sigma2 = as<double>(S_H_E_val[4]);
  
  for(int j=0; j<repm_n; j++)
  {
    m_temp = repm_times(j);
    hat_R = hat_R_full.submat(0, 0, m_temp-1, m_temp-1);  //Trim full hat_R to required m
    hee = (hat_R + (1e-6) * eye(m_temp, m_temp)).i();
    for(int k=0; k<repm(j); k++)
    {
      X = x_mat.rows(start(i), end(i));
      Y = y_vect.subvec(start(i), end(i));
      EPS = Y - X * beta_new;
      M = M + pow(sigma2, -2) * ( X.t() * hee * (EPS * EPS.t()) * hee * X );
      i++;
    }
  }
  
  hnee = (H + n * E + (1e-6) * eye(p,p)).i();
  var_sand = diagvec(hnee * M * hnee);
  
  return List::create(Named("beta.shrink") = beta_new,
                      Named("var.sand") = var_sand,
                      Named("flag.stop") = flag_stop,
                      Named("iter.n") = step_up,
                      Named("est_sigma2") = sigma2);
}

//-----------------------------------------------------
// Cross-validation across an array of lambda and pick up the best.
List CV_lam_grid_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R_full, vec beta_ini, int fold, int n, vec m, int obs_n, int p, uvec start, uvec end, vec lam_vect, double eps_stop, int iter_try){
  
  int lam_length = lam_vect.n_elem;
  double lam_temp, cv_sum, flag_stop_sum, iter_n_sum, cv_min = math::inf(), lam_min = -1;
  uvec cvgrps_seq = linspace<uvec>(0, (n-1), n); 
  uvec cvgrps_subsets = shuffle(cvgrps_seq);
  uvec cvgrps_which = cvgrps_seq - floor( cvgrps_seq / fold) * fold;
  uvec index_cv_train, index_cv_test; 
  vec cv_vect(lam_length), flag_stop_vect(lam_length), iter_n_vect(lam_length);
  uvec idx_train, idx_test;
  vec y_train, y_test, id_train, m_train, beta_train;
  mat x_train, x_test;
  List indGen_res, beta_shrink_res;
  
  for(int lam_iter = 0; lam_iter < lam_length; lam_iter++)
  {
    lam_temp = lam_vect(lam_iter);
    cv_sum = 0; 
    flag_stop_sum = 0; 
    iter_n_sum = 0;
    
    for(int k = 0; k < fold; k++)
    {
      index_cv_train = cvgrps_subsets.elem(find(cvgrps_which != k));
      index_cv_test = cvgrps_subsets.elem(find(cvgrps_which == k));
      
      idx_train = unique(seqJoin_vec(start.elem(index_cv_train), end.elem(index_cv_train), m.elem(index_cv_train)));
      idx_test = unique(seqJoin_vec(start.elem(index_cv_test), end.elem(index_cv_test), m.elem(index_cv_test)));
      
      y_train = y_vect.elem(idx_train);
      x_train = x_mat.rows(idx_train);
      id_train = id_vect.elem(idx_train);
      indGen_res = indGen_cpp(id_train);
      
      y_test = y_vect.elem(idx_test);
      x_test = x_mat.rows(idx_test);
      
      beta_shrink_res = beta_shrink_normal_cpp(y_train, x_train, id_train, hat_R_full, beta_ini, as<int>(indGen_res[0]), as<vec>(indGen_res[1]), as<int>(indGen_res[2]), p, as<uvec>(indGen_res[3]), as<uvec>(indGen_res[4]), as<vec>(indGen_res[5]), as<uvec>(indGen_res[6]), lam_temp, eps_stop, iter_try);
      cv_sum += sum(pow((y_test - x_test * as<vec>(beta_shrink_res[0])),2));                        
      flag_stop_sum += as<double>(beta_shrink_res[2]);
      iter_n_sum += as<double>(beta_shrink_res[3]);
    }
    
    // Calculate average across k-fold validation
    cv_vect(lam_iter) = cv_sum / fold;
    flag_stop_vect(lam_iter) = flag_stop_sum / fold;
    iter_n_vect(lam_iter) = iter_n_sum / fold;
    
    if(cv_sum < cv_min)
    {
      lam_min = lam_temp;
      cv_min = cv_sum;
    }
  }
  
  return List::create(Named("lam.vect") = lam_vect,
                      Named("cv.vect") = cv_vect,
                      Named("flag_stop_vect") = flag_stop_vect,
                      Named("iter_n_vect") = iter_n_vect,
                      Named("lam.min") = lam_min,
                      Named("cv.min") = cv_min
  );
}

//-----------------------------------------------------
// Best results (lambda) given Pm
// [[Rcpp::export]]
List est_pgee_grid_cpp(vec y_vect, mat x_mat, vec id_vect, vec beta_hat_R, int fold, int p, vec lam_vect, double eps_stop, int iter_try, std::string corr_str){
  
  // Initialize beta: all zeros. For CV and best model fitting, use initialized beta as starting point.
  vec beta_ini(p); beta_ini.fill(0);
  
  // Prepare sorted and unsorted index indicating repeated measures blocks
  List indGen_res = indGen_cpp(id_vect);
  
  // Given the input beta, estimate hat_R
  mat hat_R_full = corr_est_normal_cpp(y_vect, x_mat, id_vect, beta_hat_R, p, corr_str);
  
  // Run cross-validation to find the best lambda (lam.min)  
  List lam_cv = CV_lam_grid_cpp(y_vect, x_mat, id_vect, hat_R_full, beta_ini, fold, as<int>(indGen_res[0]), as<vec>(indGen_res[1]), as<int>(indGen_res[2]), p, as<uvec>(indGen_res[7]), as<uvec>(indGen_res[8]), lam_vect, eps_stop, iter_try);
  
  // Run the model using all data with the best lambda
  List beta_fit = beta_shrink_normal_cpp(y_vect, x_mat, id_vect, hat_R_full, beta_ini, as<int>(indGen_res[0]), as<vec>(indGen_res[1]), as<int>(indGen_res[2]), p, as<uvec>(indGen_res[3]), as<uvec>(indGen_res[4]), as<vec>(indGen_res[5]), as<uvec>(indGen_res[6]), as<double>(lam_cv[4]), eps_stop, iter_try);
  
  return List::create(
                      Named("lam.cv.vector") = as<vec>(lam_cv[1]),
                      Named("flag.stop.cv.vector") = as<vec>(lam_cv[2]),
                      Named("iter.cv.vector") = as<vec>(lam_cv[3]),
                      Named("lam.sel") = as<double>(lam_cv[4]),
                      Named("beta.shrink") = as<vec>(beta_fit[0]),
                      Named("var.sand") = as<vec>(beta_fit[1]),
                      Named("flag.stop") = as<double>(beta_fit[2]),
                      Named("iter.n") = as<double>(beta_fit[3]),
                      Named("est.sigma2") = as<double>(beta_fit[4]),
                      Named("hat.R") = hat_R_full
                     );
}

//-----------------------------------------------------
// Fit beta with independent structure then refit with hat_R
// [[Rcpp::export]]
List one_run_grid_cpp(vec y_vect, mat x_mat, vec id_vect, int fold, int p, vec lam_vect, double eps_stop, int iter_try, std::string corr_str){

  vec beta_ini(p); beta_ini.fill(0);
  vec beta_shrink_indep;
  
  // Find shrinked beta using independent structure, starting from initialized beta, to get shrinked beta (independent)
  List beta_fit_indep = est_pgee_grid_cpp(y_vect, x_mat, id_vect, beta_ini, fold, p, lam_vect, eps_stop, iter_try, "indep");
  beta_shrink_indep = as<vec>(beta_fit_indep[4]);
  
  // Using independent shrinked beta, estimate specific working correlation structure, and refit the model starting from initialized beta
  List beta_fit_corr = est_pgee_grid_cpp(y_vect, x_mat, id_vect, beta_shrink_indep, fold, p, lam_vect, eps_stop, iter_try, corr_str);
  
  return List::create(
                      Named("lam.cv.cor") = as<vec>(beta_fit_corr[0]),
                      Named("flag.stop.cv.cor") = as<vec>(beta_fit_corr[1]),
                      Named("iter.cv.cor") = as<vec>(beta_fit_corr[2]),
                      Named("lam.sel.cor") = as<double>(beta_fit_corr[3]),
                      Named("beta.shrink.cor") = as<vec>(beta_fit_corr[4]),
                      Named("var.sand.cor") = as<vec>(beta_fit_corr[5]),
                      Named("flag.stop.cor") = as<double>(beta_fit_corr[6]),
                      Named("iter.n.cor") = as<double>(beta_fit_corr[7]),
                      Named("est.sigma2.cor") = as<double>(beta_fit_corr[8]),
                      Named("hat.R") = as<mat>(beta_fit_corr[9])
                    );
}

