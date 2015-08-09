#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
vec q_scad_cpp(vec theta, double lambda, double a) {
  int p = theta.n_elem;
  vec b1(p), b2(p);
  b1.fill(0); b2.fill(0);
  theta = abs(theta);
  b1.elem(find(theta>lambda)).fill(1);
  b2.elem(find(theta<(lambda * a))).fill(1);
  return lambda * (1 - b1) + ((lambda * a) - theta) % b2 % (1 / (a-1) * b1);
}

// [[Rcpp::export]]
List S_H_E_normal_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R, vec beta_val, int n, int m, double lam, double eps){
  int p = x_mat.n_cols;
  double sigma = 0;
  vec beta_val_abs = abs(beta_val), Y;
  mat resi, S(p,1), H(p,p), E, X, Xh;
  S.fill(0.0); H.fill(0.0);
  
  resi = y_vect - x_mat * beta_val;
  sigma = as<double>(wrap(sqrt(sum(pow(resi, 2))/(n*m - p))));
  
  for(int i=0; i<n; i++)
  {
    X = x_mat.rows((i*m),(i*m+1));
    Y = y_vect.subvec((i*m),(i*m+1));
    Xh = X.t() * hat_R.i();
    S = S + pow(sigma,-2) * (Xh * (Y - X * beta_val));
    H = H + pow(sigma,-2) * (Xh * X);
  }
  
  E = diagmat(q_scad_cpp(beta_val_abs, lam, 3.7) / (beta_val_abs + eps));
  
  return List::create(Named("S") = S,
                      Named("H") = H,
                      Named("E") = E,
                      Named("hat.R") = hat_R,
                      Named("sigma") = sigma);
}

// [[Rcpp::export]]
List beta_shrink_normal_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R, vec beta_ini, int n, int m, double lam, double eps, int eps_stop, int max_step){
  int p = x_mat.n_cols;
  vec beta_new = beta_ini, Y, var_sand;
  int step_up = 0;
  mat S, H, E, M(p,p), X, EPS, hee, hnee;
  M.fill(0.0);
  double flag_stop, sigma;
  
  List S_H_E_val = S_H_E_normal_cpp(y_vect, x_mat, id_vect, hat_R, beta_new, n, m, lam, eps);
  S = as<mat>(S_H_E_val[0]);
  H = as<mat>(S_H_E_val[1]);
  E = as<mat>(S_H_E_val[2]);
  while(step_up < max_step)
  {
    vec beta_old = beta_new;
    beta_new = beta_old + (H + n * E + (1e-6) * eye(p,p)).i() * (S - n * E * beta_old);
    S_H_E_val = S_H_E_normal_cpp(y_vect, x_mat, id_vect, hat_R, beta_new, n, m, lam, eps);
    S = as<mat>(S_H_E_val[0]);
    H = as<mat>(S_H_E_val[1]);
    E = as<mat>(S_H_E_val[2]);
    step_up++;
    flag_stop = sum(abs(beta_old-beta_new));
  }
  
  hat_R = as<mat>(S_H_E_val[3]);
  sigma = as<double>(S_H_E_val[4]);
  hee = (hat_R + (1e-6) * eye(m,m)).i();
  
  for(int i=0; i<n; i++)
  {
    X = x_mat.rows((i*m),(i*m+1));
    Y = y_vect.subvec((i*m),(i*m+1));
    EPS = Y - X * beta_new;
    M = M + pow(sigma, -4) * ( X.t() * hee * (EPS * EPS.t()) * hee * X );
  }
  
  hnee = (H + n * E + (1e-6) * eye(p,p)).i();
  
  var_sand = diagvec(hnee * M * hnee);
  
  return List::create(Named("beta.shrink") = beta_new,
                      Named("beta.ini") = beta_ini,
                      Named("var.sand") = var_sand,
                      Named("flag.stop") = flag_stop,
                      Named("est_sigma2") = pow(sigma,2));
}


// [[Rcpp::export]]
List CV_lam_grid_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R, vec beta_ini, int n, int m, int fold, vec lam_vect, double eps, int eps_stop, int max_step){
  int lam_length = lam_vect.n_elem;
  int n_train;
  double lam_temp, cv_value, cv_min = math::inf(), lam_min = -1;
  vec cvgrps_seq = linspace<vec>(0, (n-1), n); 
  vec cvgrps_seqm = m*cvgrps_seq;
  vec cvgrps_subsets = shuffle(cvgrps_seqm);
  vec cvgrps_which = cvgrps_seq - floor( cvgrps_seq / fold) * fold;
  vec index_cv_train, index_cv_1_train, index_cv_test, index_cv_1_test, cv_vect(lam_length); cv_vect.fill(0);
  uvec idx_train, idx_test;
  vec y_train, y_test, id_train, id_test, beta_train;
  mat x_train, x_test;
  
  for(int j=0; j<lam_length; j++)
  {
      lam_temp = lam_vect(j);
      cv_value = 0;
      for(int k=0; k<fold; k++)
      {
          index_cv_train = index_cv_1_train =cvgrps_subsets.elem(find(cvgrps_which!=k));
          index_cv_test = index_cv_1_test =cvgrps_subsets.elem(find(cvgrps_which==k));
  
            for (int i=0; i<(m-1); i++)
            {
              idx_train = as<uvec>(wrap(join_cols((index_cv_1_train + (i + 1)) , index_cv_train)));
              idx_test = as<uvec>(wrap(join_cols((index_cv_1_test + (i + 1)) , index_cv_test)));
            }
          y_train = y_vect.elem(idx_train);
          x_train = x_mat.rows(idx_train);
          id_train = id_vect.elem(idx_train);
          n_train = idx_train.n_elem / m;
            
          y_test = y_vect.elem(idx_test);
          x_test = x_mat.rows(idx_test);
          id_test = id_vect.elem(idx_test);
          
          beta_train = as<vec>(beta_shrink_normal_cpp(y_train, x_train, id_train, hat_R, beta_ini, n_train, m, lam_temp, eps, eps_stop, max_step)[0]);
          cv_value += sum(pow((y_test - x_test * beta_train),2));                        
      }
      cv_vect(j) = cv_value;
      if(cv_value < cv_min)
      {
        lam_min = lam_temp;
        cv_min = cv_value;
      }
  }
  return List::create(Named("lam.vect") = lam_vect,
                      Named("cv.vect") = cv_vect,
                      Named("lam.min") = lam_min,
                      Named("cv.min") = cv_min);
}

// [[Rcpp::export]]
List est_pgee_grid_cpp(vec y_vect, mat x_mat, vec id_vect, mat hat_R, vec beta_ini, int n, int m, int fold, vec lam_vect, double eps, int eps_stop, int max_step){
  List lam_cv = CV_lam_grid_cpp(y_vect, x_mat, id_vect, hat_R, beta_ini, n, m, fold, lam_vect, eps, eps_stop, max_step);
  double lam_sel = as<double>(lam_cv[2]);
  // Run the model using all data with the best lambda
  List beta_fit = beta_shrink_normal_cpp(y_vect, x_mat, id_vect, hat_R, beta_ini, n, m, lam_sel, eps, eps_stop, max_step);
  return List::create(
                      Named("lam.cv.vector") = as<vec>(lam_cv[1]),
                      Named("beta.shrink") = as<vec>(beta_fit[0]),
                      Named("var.sand") = as<vec>(beta_fit[2]),
                      Named("flag.stop") = as<double>(beta_fit[3]),
                      Named("lam.sel") = lam_sel,
                      Named("est.sigma2") = as<double>(beta_fit[4])
                      );
}

// [[Rcpp::export]]

mat corr_est_normal_cpp(vec y_vect, mat x_mat, vec id_vect, vec beta_val, int n, int m, std::string corr_str){
  int p = x_mat.n_cols;
  vec resi = y_vect - x_mat * beta_val, temp_resi, ones(m), alphas(m); ones.fill(1);
  double sigma = as<double>(wrap(sqrt(sum(pow(resi, 2))/(n * m - p))));
  mat hat_R(m,m), resi_mat, hat_l(m,m), hat_r(m,m); hat_R.fill(0); hat_r.fill(0);
  double SUM, SUM_i, alpha;
  resi = resi / sigma;
  
  if(corr_str == "indep")
  {
    hat_R = eye(m,m);
  }
  
  if(corr_str == "exch")
  {
    SUM = 0;
    for (int i=0; i<n; i++)
    {
      temp_resi = resi.subvec((i*m),(i*m+1));
      SUM_i = 0;
      for (int j=0; j<m; j++)
      {
        SUM_i = SUM_i + (sum(temp_resi(j)*temp_resi) - pow(temp_resi(j),2));
      }
      SUM = SUM + SUM_i / (m * (m-1));
    }
    alpha = SUM / n;
    alphas.fill(1-alpha);
    hat_r.diag() = alphas;
    hat_R = hat_l.fill(alpha) + hat_r;
  }
  
  if (corr_str == "ar1")
  {
    SUM = 0;
    for (int i=0; i<n; i++)
    {
      temp_resi = resi.subvec((i*m), (i*m+1));
      SUM = SUM + sum(temp_resi.subvec(0, (m-2)) * temp_resi.subvec(1, (m-1))) / (m - 1);
    }
    alpha = SUM / n;
    for (int j=0; j<m; j++ )
    {
      for (int k=0; k<m; k++)
      {
        hat_R(j,k) = pow(alpha,abs((k + 1) - (j + 1)));
      }
    }
  }
  
  if (corr_str == "un")
  {
    resi_mat.insert_cols(0, resi);
    resi_mat.reshape(m, n);
    resi_mat = resi_mat.t();
    hat_R = resi_mat.t() * resi_mat / n;
    hat_R.diag() = ones;
  }
  return(hat_R);
}


// [[Rcpp::export]]
List one_run_grid_cpp(vec y_vect, mat x_mat, vec id_vect, int n, int m, int pn, int fold, vec lam_vect, double rho, double eps, int eps_stop, int max_step, std::string corr_str){
  mat hat_ini(m,m), hat_R_corr; hat_ini.eye();
  vec beta_ini(pn); beta_ini.fill(0);
  vec beta_shrink;
  List beta_fit_indep = est_pgee_grid_cpp(y_vect, x_mat, id_vect, hat_ini, beta_ini, n, m, fold, lam_vect, eps, eps_stop, max_step);
  beta_shrink = as<vec>(beta_fit_indep[1]);
  //Estimate the correlation structure within subject and refit the model
  hat_R_corr = corr_est_normal_cpp(y_vect, x_mat, id_vect, beta_shrink, n, m, corr_str);
  List beta_fit_corr = est_pgee_grid_cpp(y_vect, x_mat, id_vect, hat_R_corr, beta_ini, n, m, fold, lam_vect, eps, eps_stop, max_step);
  return List::create(
    Named("lam.cv.vector") = as<vec>(beta_fit_corr[0]),
    Named("beta.shrink.corr") = as<vec>(beta_fit_corr[1]),
    Named("var.sand.corr") = as<vec>(beta_fit_corr[2]),
    Named("flag.stop.corr") = as<double>(beta_fit_corr[3]),
    Named("lam.sel.corr") = as<double>(beta_fit_corr[4]),
    Named("est.sigma2") = as<double>(beta_fit_corr[5]),
    Named("alpha.corr") = hat_R_corr
  );
}