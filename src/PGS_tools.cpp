#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//-----------------------------------------------------
// Count frequencies of array element
// [[Rcpp::export]]
vec countRep_cpp(vec id){
  
  int id_len = id.n_elem;
  vec buffer(1); buffer.fill(-1);
  id = join_cols(id, buffer);
  vec count_vec(id_len); count_vec.fill(0);
  int prev = id(0);
  int count = 1, j = 0;
  for (int i = 1; i <= id_len; i++) {
    if (id(i) == prev) {
      count++;
    } else {
      prev = id(i);
      count_vec(j) = count;
      count = 1;
      j++;
    }
  }
  count_vec = count_vec.elem(find(count_vec>0));
  return(count_vec);
}

//-----------------------------------------------------
// Generate index block for same number of repeated measures
// [[Rcpp::export]]
List indGen_cpp(vec id_vect){
  
  vec m = countRep_cpp(id_vect);
  vec end = cumsum(m) - 1;
  vec start = end - m + 1;
  uvec sort_ind = sort_index(m);
  vec start_sort = start.elem(sort_ind);
  vec end_sort = end.elem(sort_ind);
  
  return List::create(Named("n") = m.n_elem,
                      Named("m") = m,
                      Named("obs_n") = sum(m),
                      Named("start_sort") = start_sort,
                      Named("end_sort") = end_sort,
                      Named("repm") = countRep_cpp(m(sort_ind)),
                      Named("repm_times") = unique(m),
                      Named("start") = start,
                      Named("end") = end);
}

//-----------------------------------------------------
// Join sequence of elements of two vector:  m is a vector (for unsorted seq)
// [[Rcpp::export]]
uvec seqJoin_vec(uvec seq1, uvec seq2, vec m){
  
  int seq_len = seq1.n_elem;
  uvec res(sum(m)); res.fill(0);
  uvec temp;
  int j = 0;
  
  for (int i = 0; i < seq_len; i++)
  {
    temp = linspace<uvec>(seq1(i),seq2(i), m(i));
    res.subvec(j, j + temp.n_elem - 1) = temp;
    j += temp.n_elem;
  }
  
  return(res);
}

//-----------------------------------------------------
// Join sequence of elements of two vector: m is an integer (for sorted seq)
// [[Rcpp::export]]
uvec seqJoin_int(uvec seq, int m){
  
  uvec seq1 = seq;
  for (int i = 1; i < m; i++)
  {
    seq = join_cols(seq, seq1 + i);
  }
  seq = sort(seq);
  
  return(seq);
}

//-----------------------------------------------------
// Estimate working correlation matrix
// [[Rcpp::export]]
mat corr_est_normal_cpp(vec y_vect, mat x_mat, vec id_vect, vec beta_val, int p, std::string corr_str){
  
  //Extract sorted index information 
  List indGen_res = indGen_cpp(id_vect);
    int n = as<int>(indGen_res[0]);
    vec m = as<vec>(indGen_res[1]);
    int obs_n = as<int>(indGen_res[2]);
    uvec start = as<uvec>(indGen_res[3]);
    uvec end = as<uvec>(indGen_res[4]);
    vec repm = as<vec>(indGen_res[5]);
    uvec repm_times = as<uvec>(indGen_res[6]);
  
  int i = 0, m_max = max(m), repm_n = repm.n_elem, m_temp;
  vec temp_resi, ones(m_max), alphas(m_max); 
  double SUM, SUM_i, alpha;
  mat hat_R(m_max, m_max), hat_l(m_max, m_max), hat_r(m_max, m_max); hat_R.fill(0); 
  uvec idx_repm, repm_seq;
  vec repm_cum, buffer(1), div(m_max);
  
  vec resi = y_vect - x_mat * beta_val;
  resi =  resi / as<double>(wrap(sqrt(sum(pow(resi, 2))/(obs_n - p))));
  
  // Independent
  if(corr_str == "indep")
  {
    hat_R = eye(m_max, m_max);
  }
  
  // Exchangable
  if(corr_str == "exch")
  {
    hat_r.fill(0);
    SUM = 0;
    for (int j = 0; j < repm_n; j++)
    {
      m_temp = repm_times(j);
      for (int k = 0; k < repm(j); k++)
      {
        temp_resi = resi.subvec(start(i), end(i));
        SUM_i = 0;
        if (m_temp > 1)
        {
          for (int repi = 0; repi < m_temp; repi++)
          {
            SUM_i = SUM_i + (sum(temp_resi(repi)*temp_resi) - pow(temp_resi(repi),2));
          }
          SUM += SUM_i / (m_temp * (m_temp - 1));
        }
        i++;
      }
    }
    alpha = SUM / n;
    alphas.fill(1-alpha);
    hat_r.diag() = alphas;
    hat_R = hat_l.fill(alpha) + hat_r;
  }
  
  // Auto-regressive(1)
  if (corr_str == "ar1")
  {
    SUM = 0;
    for (int j = 0; j < repm_n; j++)
    {
      m_temp = repm_times(j);
      for (int k = 0; k < repm(j); k++)
      {
        if (m_temp > 1)
        {
          temp_resi = resi.subvec(start(i), end(i));
          SUM += sum(temp_resi.subvec(0, (m_temp - 2)) % temp_resi.subvec(1, (m_temp - 1))) / (m_temp - 1);
        }
        i++;
      }
    }
    alpha = SUM / n;
    
    for (int elemj = 0; elemj < m_max; elemj++ )
    {
      for (int elemk = 0; elemk < m_max; elemk++)
      {
        hat_R(elemj, elemk) = pow(alpha, abs((elemk + 1) - (elemj + 1)));
      }
    }
  }
  
  //Unstructured
  if (corr_str == "un")
  {
    buffer.fill(0);
    div.fill(0);
    ones.fill(1);
    repm_cum = join_cols(buffer, cumsum(repm)) - 1;
    repm_n = repm_cum.n_elem;
    for (int j = 0; j < repm_n - 1; j++)
    {
      mat resi_mat, hat_R_temp;
      m_temp = repm_times(j);
      if(m_temp > 1)
      {
        repm_seq = linspace<uvec>(repm_cum(j)+1, repm_cum(j+1), repm(j)); 
        idx_repm = seqJoin_int(start.elem(repm_seq), m_temp);
        temp_resi = resi.elem(idx_repm);
        resi_mat.insert_cols(0, temp_resi);
        resi_mat.reshape(m_temp, repm(j));
        hat_R_temp = resi_mat * resi_mat.t();
        hat_R_temp.resize(m_max, m_max);
        hat_R += hat_R_temp;
      }
    }
    
    div.elem(repm_times-1) = repm;
    std::reverse(div.begin(), div.end());
    div = cumsum(div);
    std::reverse(div.begin(), div.end());
    
    for (int hati = 1; hati < m_max; hati++)
    {
      for(int hatj = 0; hatj < hati; hatj++)
      {
        hat_R(hati, hatj) = hat_R(hati, hatj) / div(hati);
        hat_R(hatj, hati) = hat_R(hatj, hati) / div(hati);          
      }
    }
    
    hat_R.diag() = ones;
  }
  
  //   return List::create(
  //     Named("repm_cum") = repm_cum,
  //     Named("repm_n") = repm_n,
  //     Named("repm_seq") = repm_seq,
  //     Named("div") = div,
  //     Named("hat_R") = hat_R
  //   );
  return(hat_R);
}
