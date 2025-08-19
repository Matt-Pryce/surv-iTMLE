// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <chrono>
using namespace Rcpp;
using namespace arma;


ivec findInterval_cpp(const vec& x, const vec& v) {
  int nx = x.n_elem, nv = v.n_elem;
  ivec res(nx);
  for (int i = 0; i < nx; ++i) {
    int j = 0;
    while (j < nv && x[i] >= v[j]) ++j;
    res[i] = j;
  }
  return res;
}


// [[Rcpp::export]]
mat int_fmx_dF_cpp(const vec& v, const mat& fmx, const mat& Fmx, const vec& jumps) {
  int n = fmx.n_rows;
  int m = fmx.n_cols;
  
  // Step 1: Compute discrete delta of CDF
  mat dFmx = Fmx;
  dFmx.cols(1, m-1) -= Fmx.cols(0, m-2);
  dFmx.col(0) = Fmx.col(0);  // leftmost jump
  
  // Step 2: Elementwise product
  mat vals = fmx % dFmx;
  
  // Step 3: Prefix sum across columns (row-wise)
  mat cumsum_vals = cumsum(vals, 1);
  
  // Step 4: Lookup integral by jump index
  ivec id = findInterval_cpp(v, jumps);
  mat result(n, v.n_elem, fill::zeros);
  
  for (int k = 0; k < v.n_elem; ++k) {
    if (id[k] == 0) continue;  // integral is 0
    result.col(k) = cumsum_vals.col(id[k] - 1);
  }
  
  return result;
}


// // [[Rcpp::export]]
// mat int_infty_fmx_dF_cpp(const vec& v, const mat& fmx, const mat& Fmx, const vec& jumps) {
//   int n = fmx.n_rows;
//   int m = fmx.n_cols;
//   mat dFmx = Fmx;
//   dFmx.cols(1, m-1) -= Fmx.cols(0, m-2);
//   dFmx.col(0) = Fmx.col(0);
//   mat vals = fmx % dFmx;
//   mat result(n, v.n_elem, fill::zeros);
//   ivec id = findInterval_cpp(v, jumps);
// 
//   for (int k = 0; k < v.n_elem; ++k) {
//     for (int i = 0; i < n; ++i) {
//       for (int j = id[k]; j < m; ++j) {
//         result(i, k) += vals(i, j);
//       }
//     }
//   }
//   return result;
// }


// [[Rcpp::export]]
mat int_infty_fmx_dF_cpp(const vec& v, const mat& fmx, const mat& Fmx, const vec& jumps) {
  int n = fmx.n_rows;
  int m = fmx.n_cols;
  
  // Step 1: Compute dF
  mat dFmx = Fmx;
  dFmx.cols(1, m - 1) -= Fmx.cols(0, m - 2);
  dFmx.col(0) = Fmx.col(0);
  
  // Step 2: Weighted integrand
  mat vals = fmx % dFmx;
  
  // Step 3: Reverse columns
  umat rev_indices = linspace<umat>(m - 1, 0, m);
  mat vals_rev = vals.cols(rev_indices);
  
  // Step 4: Row-wise cumsum on reversed matrix
  mat cumsum_rev = cumsum(vals_rev, 1);
  
  // Step 5: Flip back to get reverse cumsum
  mat rev_cumsum_vals = cumsum_rev.cols(rev_indices);
  
  // Step 6: Fast lookup
  mat result(n, v.n_elem, fill::zeros);
  ivec id = findInterval_cpp(v, jumps);
  
  for (int k = 0; k < v.n_elem; ++k) {
    if (id[k] < m) {
      result.col(k) = rev_cumsum_vals.col(id[k]);
    }
  }
  
  return result;
}




// [[Rcpp::export]]
mat CDF_eval_mx_cpp(const vec& time_eval, const mat& CDF_mx, const vec& jumps) {
  int n = CDF_mx.n_rows;
  mat CDF_mx_ext(n, CDF_mx.n_cols + 1, fill::zeros);
  CDF_mx_ext.cols(1, CDF_mx.n_cols) = CDF_mx;
  mat out(n, time_eval.n_elem);

  ivec id = findInterval_cpp(time_eval, jumps);
  for (int k = 0; k < time_eval.n_elem; ++k) {
    out.col(k) = CDF_mx_ext.col(id[k]);
  }
  return out;
}

// // Linear search algorithm
// mat CDF_eval_mx_cpp_matrix(const mat& time_eval_mx, const mat& CDF_mx, const vec& jumps) {
//   int n = CDF_mx.n_rows;
//   int m = CDF_mx.n_cols;
//   
//   mat CDF_mx_ext(n, m + 1, fill::zeros);
//   CDF_mx_ext.cols(1, m) = CDF_mx;
//   
//   mat out(n, m);
//   
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < m; ++j) {
//       double t = time_eval_mx(i, j);
//       int id = 0;
//       while (id < jumps.n_elem && t >= jumps[id]) ++id;
//       out(i, j) = CDF_mx_ext(i, id);
//     }
//   }
//   return out;
// }

// // Helper: Binary search to find interval
// inline int find_interval_binary(double t, const vec& jumps) {
//   int left = 0, right = jumps.n_elem;
//   while (left < right) {
//     int mid = (left + right) / 2;
//     if (t >= jumps[mid]) {
//       left = mid + 1;
//     } else {
//       right = mid;
//     }
//   }
//   return left;  // returns index where t < jumps[id]
// }

// // Using binary search, which is faster
// mat CDF_eval_mx_cpp_matrix(const mat& time_eval_mx, const mat& CDF_mx, const vec& jumps) {
//   
//   int n = CDF_mx.n_rows;
//   int m = CDF_mx.n_cols;
//   
//   // Extend CDF matrix with a 0 column at the beginning
//   mat CDF_mx_ext(n, m + 1, fill::zeros);
//   CDF_mx_ext.cols(1, m) = CDF_mx;
//   
//   // Output matrix
//   mat out(n, time_eval_mx.n_cols);
//   
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < time_eval_mx.n_cols; ++j) {
//       double t = time_eval_mx(i, j);
//       int id = find_interval_binary(t, jumps);   // Using binary search, which is faster
//       out(i, j) = CDF_mx_ext(i, id);  // id is already 0-indexed
//     }
//   }
//   
//   return out;
// }


mat CDF_eval_mx_cpp_matrix_rowIncreasing(const mat& time_eval_mx, const mat& CDF_mx, const vec& jumps) {
  // Suppose mat is non-decreasing for each row.
  int n = CDF_mx.n_rows;
  int m = CDF_mx.n_cols;
  int k = jumps.n_elem;
  
  // Extend CDF matrix with a 0 column at the beginning
  mat CDF_mx_ext(n, m + 1, fill::zeros);
  CDF_mx_ext.cols(1, m) = CDF_mx;
  
  // Output matrix
  mat out(n, time_eval_mx.n_cols);
  
  for (int i = 0; i < n; ++i) {
    int id = 0;
    for (int j = 0; j < time_eval_mx.n_cols; ++j) {
      double t = time_eval_mx(i, j);
      while (id < k && t >= jumps[id]) ++id;
      out(i, j) = CDF_mx_ext(i, id);
    }
  }
  
  return out;
}



// mat int_infty_fmx_dF_batch_cpp(const mat& v_mat,
//                                const cube& f_cube,
//                                const cube& F_cube,
//                                const vec& jumps) {
//   int n = v_mat.n_rows;
//   int m = v_mat.n_cols;
//   int l = jumps.n_elem;
//   
//   mat result(n, m, fill::zeros);
//   
//   // Compute dF, the 3rd dimension (slice) is the time dimension
//   cube dF = F_cube;
//   for (int i = 0; i < F_cube.n_slices; ++i) {
//     if (i == 0) {
//       dF.slice(i) = F_cube.slice(i);
//     } else {
//       dF.slice(i) = F_cube.slice(i) - F_cube.slice(i - 1);
//     }
//   }
//   
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < m; ++j) {
//       double vij = v_mat(i, j);
//       int idx = 0;
//       while (idx < l && vij >= jumps[idx]) ++idx;
//       
//       for (int k = idx; k < l; ++k) {
//         result(i, j) += f_cube(i, j, k) * dF(i, j, k);
//       }
//     }
//   }
//   
//   return result;
// }





// [[Rcpp::export]]
List aug_QD_1_cpp1(int nn,
                   const vec& nuu,
                   const vec& X,
                   const vec& Q,
                   const ivec& Delta,
                   const mat& Fuz_mx,
                   const mat& Gvz_mx,
                   const mat& Sdz_mx,
                   const vec& u,
                   const vec& v,
                   const vec& d,
                   double trim = 1e-7) {
  
  int u_len = u.n_elem;
  int d_len = d.n_elem;
  
  vec X_res = X - Q;
  double tau_Tmax = u.max() + 1e-5;
  double tau_Dmax = d.max() + 1e-5;
  
  // IPCW weights
  mat one_minus_Sdz = 1.0 - Sdz_mx;
  vec Sdyz_vec = 1.0 - CDF_eval_mx_cpp(X_res, one_minus_Sdz, d).diag();
  vec w_D = zeros<vec>(nn);
  uvec id0 = find(Delta == 0);
  w_D(id0) = 1.0 / clamp(Sdyz_vec(id0), trim, datum::inf);
  
  // aug_11 and aug_11_const1
  vec Fxz_vec = CDF_eval_mx_cpp(X, Fuz_mx, u).diag();
  mat Guz_mx = CDF_eval_mx_cpp(u, Gvz_mx, v);  // nn x u
  
  mat nuu_mx = repmat(nuu.t(), nn, 1);
  mat ind_Xu_mx = conv_to<mat>::from(repmat(X, 1, u_len) <= repmat(u.t(), nn, 1));
  mat Guz_trim = clamp(Guz_mx, trim, datum::inf);
  
  mat fu1_mx = nuu_mx % ind_Xu_mx / Guz_trim;
  vec aug_11 = w_D % int_fmx_dF_cpp(vec({tau_Tmax}), fu1_mx, Fuz_mx, u).col(0) /
    clamp(1.0 - Fxz_vec, trim, datum::inf);
  
  mat fu1_const1_mx = ind_Xu_mx / Guz_trim;
  vec aug_11_const1 = w_D % int_fmx_dF_cpp(vec({tau_Tmax}), fu1_const1_mx, Fuz_mx, u).col(0) /
    clamp(1.0 - Fxz_vec, trim, datum::inf);
  
  // aug_12 and aug_12_const1
  mat wnu_mx = nuu_mx / Guz_trim;
  mat wnu_const1_mx = 1.0 / Guz_trim;
  mat int_wnu_dF(nn, d_len, fill::zeros);
  mat int_wnu_const1_dF(nn, d_len, fill::zeros);
  
  mat Q_mat = repmat(Q, 1, d_len);
  mat d_mat = repmat(d.t(), nn, 1);
  mat Q_plus_d_mx = Q_mat + d_mat;
  mat F_Qu_z_mx = CDF_eval_mx_cpp_matrix_rowIncreasing(Q_plus_d_mx, Fuz_mx, u);
  mat F_Qu_z_trim = clamp(1.0 - F_Qu_z_mx, trim, datum::inf);
  
  vec Q_plus_d = vec(d_len);  // reuse buffer
  
  for (int i = 0; i < nn; ++i) {
    Q_plus_d = Q[i] + d;  // inplace addition
    rowvec Fuz_row = Fuz_mx.row(i);
    rowvec wnu_row = wnu_mx.row(i);
    rowvec wnu_const1_row = wnu_const1_mx.row(i);
    
    int_wnu_dF.row(i) = int_infty_fmx_dF_cpp(Q_plus_d, wnu_row, Fuz_row, u);
    int_wnu_const1_dF.row(i) = int_infty_fmx_dF_cpp(Q_plus_d, wnu_const1_row, Fuz_row, u);
  }
  
  mat X_res_mat = repmat(X_res, 1, d_len);
  mat ind_dXQ_mx = conv_to<mat>::from(d_mat <= X_res_mat);
  
  mat Sdz_trim = clamp(Sdz_mx, trim, datum::inf);
  mat FDdz_mx = 1.0 - Sdz_mx;
  mat Sdz_trim_sq = Sdz_trim % Sdz_trim;
  
  mat fd2_mx = ind_dXQ_mx % int_wnu_dF / (F_Qu_z_trim % Sdz_trim_sq);
  vec aug_12 = int_fmx_dF_cpp(vec({tau_Dmax}), fd2_mx, FDdz_mx, d).col(0);
  
  mat fd2_const1_mx = ind_dXQ_mx % int_wnu_const1_dF / (F_Qu_z_trim % Sdz_trim_sq);
  vec aug_12_const1 = int_fmx_dF_cpp(vec({tau_Dmax}), fd2_const1_mx, FDdz_mx, d).col(0);
  
  vec aug_1 = aug_11 - aug_12;
  vec aug_1_const1 = aug_11_const1 - aug_12_const1;
  
  return List::create(Named("aug_1") = aug_1,
                      Named("aug_1_const1") = aug_1_const1);
}






// [[Rcpp::export]]
List aug_QD_2_cpp1(int nn,
                  const vec& nuu,
                  const vec& Q,
                  const ivec& Delta,
                  const mat& Fuz_mx,
                  const mat& Gvz_mx,
                  const vec& u,
                  const vec& v,
                  const vec& Sdyz_vec,
                  double trim) {
  
  // replicate nuu across rows
  mat nuu_mx = repmat(nuu.t(), nn, 1); // nn × length(u)
  
  // int_fmx_dF(Q, nuu_mx, Fuz_mx) returns nn × 1 matrix, take diagonal
  vec mqz = int_fmx_dF_cpp(Q, nuu_mx, Fuz_mx, u).diag();  // vector of length nn
  
  // Gqz = G(Q|A,Z), Fqz = F(Q|A,Z)
  vec Gqz = CDF_eval_mx_cpp(Q, Gvz_mx, v).diag();
  vec Fqz = CDF_eval_mx_cpp(Q, Fuz_mx, u).diag();
  
  // aug_2 = (1 - Delta / Sdyz) * mqz / [(1 - Fqz) * Gqz]
  vec denom = clamp(1.0 - Fqz, trim, datum::inf) % clamp(Gqz, trim, datum::inf);
  vec aug_2 = (1.0 - conv_to<vec>::from(Delta) / clamp(Sdyz_vec, trim, datum::inf)) % mqz / denom;
  
  // aug_2_const1 = (1 - Delta / Sdyz) * Fqz / [(1 - Fqz) * Gqz]
  vec aug_2_const1 = (1.0 - conv_to<vec>::from(Delta) / clamp(Sdyz_vec, trim, datum::inf)) % Fqz / denom;
  
  return List::create(
    Named("aug_2") = aug_2,
    Named("aug_2_const1") = aug_2_const1
  );
}



// [[Rcpp::export]]
List aug_QD_3_cpp1(int nn,
                   const vec& nuu,
                   const vec& X,
                   const vec& Q,
                   const ivec& Delta,
                   const mat& Fuz_mx,
                   const mat& Gvz_mx,
                   const mat& Sdz_mx,
                   const vec& u,
                   const vec& v,
                   const vec& d,
                   double trim) {

  int v_len = v.n_elem;
  int d_len = d.n_elem;

  vec X_res = X - Q;
  double tau_Qmax = v.max() + 1e-5;
  double tau_Dmax = d.max() + 1e-5;

  // IPCW weights
  mat one_minus_Sdz = 1.0 - Sdz_mx;
  vec Sdyz_vec = 1.0 - CDF_eval_mx_cpp(X_res, one_minus_Sdz, d).diag();
  vec w_D(nn, fill::zeros);
  uvec id0 = find(Delta == 0);
  w_D(id0) = 1.0 / clamp(Sdyz_vec(id0), trim, datum::inf);

  // aug_31
  mat nuu_mx = repmat(nuu.t(), nn, 1);
  mat mvz_mx = int_fmx_dF_cpp(v, nuu_mx, Fuz_mx, u);
  mat Fvz_mx = CDF_eval_mx_cpp(v, Fuz_mx, u);

  // Precompute elementwise min(X, v)
  mat X_mat = repmat(X, 1, v_len);
  mat v_mat = repmat(v.t(), nn, 1);
  mat min_Xv_mx = arma::min(X_mat, v_mat);
  mat F_Xv_z_mx = CDF_eval_mx_cpp_matrix_rowIncreasing(min_Xv_mx, Fuz_mx, u);

  // Precompute indicators
  mat Q_mat = repmat(Q, 1, v_len);
  mat ind_Qv_mx = conv_to<mat>::from(Q_mat <= v_mat);

  mat Gvz_trim = clamp(Gvz_mx, trim, datum::inf);
  mat Gvz_trim_sq = Gvz_trim % Gvz_trim;
  mat one_minus_F_Xv = clamp(1.0 - F_Xv_z_mx, trim, datum::inf);

  mat fv31_mx = mvz_mx % ind_Qv_mx / (one_minus_F_Xv % Gvz_trim_sq);
  vec aug_31 = w_D % int_fmx_dF_cpp(vec({tau_Qmax}), fv31_mx, Gvz_mx, v).col(0);

  mat fv31_const1_mx = Fvz_mx % ind_Qv_mx / (one_minus_F_Xv % Gvz_trim_sq);
  vec aug_31_const1 = w_D % int_fmx_dF_cpp(vec({tau_Qmax}), fv31_const1_mx, Gvz_mx, v).col(0);

  // aug_32
  mat aug32_int_dG(nn, d_len, fill::zeros);
  mat aug32_int_dG_const1(nn, d_len, fill::zeros);

  mat mvz_rep(d_len, v_len);
  mat Fvz_rep(d_len, v_len);
  mat Gvz_rep(d_len, v_len);
  mat ind_Qv_rep(d_len, v_len);

  mat Qd_mat = repmat(Q, 1, d_len);
  mat d_mat = repmat(d.t(), nn, 1);
  mat Q_plus_d_mx = Qd_mat + d_mat;

  // Precompute v_mat_d once (d_len × v_len)
  mat v_mat_d = repmat(v.t(), d_len, 1);
  
  // Preallocate reused buffers
  mat Qdv_mx(d_len, v_len);
  mat eval_grid(d_len, v_len);
  mat F_Qdv_z(d_len, v_len);
  mat fv32(d_len, v_len);
  mat fv32_const1(d_len, v_len);
  mat denom(d_len, v_len);
  
  for (int i = 0; i < nn; ++i) {
    // Construct Q + d row
    vec Q_plus_d = Q_plus_d_mx.row(i).t();
    Qdv_mx.each_col() = Q_plus_d;
    eval_grid = arma::min(Qdv_mx, v_mat_d);
    
    // auto t1 = std::chrono::high_resolution_clock::now();  // keep track of time spent
    
    // Inlined, row-wise rolling search with reset id
    {
      rowvec Fuz_row = Fuz_mx.row(i);
      rowvec Fuz_ext(u.n_elem + 1, fill::zeros);
      Fuz_ext.cols(1, u.n_elem) = Fuz_row;
      
      for (int r = 0; r < d_len; ++r) {
        int id = 0;  // 🔧 reset index for each row!
        for (int c = 0; c < v_len; ++c) {
          double t = eval_grid(r, c);
          while (id < u.n_elem && t >= u[id]) ++id;
          F_Qdv_z(r, c) = Fuz_ext[id];
        }
      }
    }
    
    // auto t2 = std::chrono::high_resolution_clock::now(); // keep track of time spent
    
    // Reuse precomputed row slices
    rowvec mvz_row = mvz_mx.row(i);
    rowvec Fvz_row = Fvz_mx.row(i);
    rowvec Gvz_row = Gvz_mx.row(i);
    
    mvz_rep.each_row() = mvz_row;
    Fvz_rep.each_row() = Fvz_row;
    Gvz_rep.each_row() = Gvz_row;
    ind_Qv_rep.each_row() = ind_Qv_mx.row(i);
    
    // Denominator computation
    denom = clamp(1.0 - F_Qdv_z, trim, datum::inf);
    denom %= Gvz_rep;
    denom %= Gvz_rep;
    
    // First integral
    fv32 = mvz_rep;
    fv32 %= ind_Qv_rep;
    fv32 /= denom;
    aug32_int_dG.row(i) = int_fmx_dF_cpp(vec({tau_Qmax}), fv32, Gvz_rep, v).col(0).t();
    
    // Second integral
    fv32_const1 = Fvz_rep;
    fv32_const1 %= ind_Qv_rep;
    fv32_const1 /= denom;
    aug32_int_dG_const1.row(i) = int_fmx_dF_cpp(vec({tau_Qmax}), fv32_const1, Gvz_rep, v).col(0).t();
    
    // auto t3 = std::chrono::high_resolution_clock::now();  // keep track of time spent
    
    // // Timing output at last iteration
    // if (i == nn - 1) {
    //   std::chrono::duration<double> loop_duration_CDFeval = t2 - t1;
    //   Rcpp::Rcout << "Time for CDF_eval: " << loop_duration_CDFeval.count() << " seconds" << std::endl;
    //   std::chrono::duration<double> loop_duration_int = t3 - t2;
    //   Rcpp::Rcout << "Time for the two int: " << loop_duration_int.count() << " seconds" << std::endl;
    // }
  }
  

  mat X_res_mat = repmat(X_res, 1, d_len);
  mat atrisk_XQ = conv_to<mat>::from(X_res_mat >= d_mat);

  mat FDdz_mx = 1.0 - Sdz_mx;
  mat Sdz_trim_sq = pow(clamp(Sdz_mx, trim, datum::inf), 2);

  mat fd32 = aug32_int_dG % atrisk_XQ / Sdz_trim_sq;
  vec aug_32 = int_fmx_dF_cpp(vec({tau_Dmax}), fd32, FDdz_mx, d).col(0);

  mat fd32_const1 = aug32_int_dG_const1 % atrisk_XQ / Sdz_trim_sq;
  vec aug_32_const1 = int_fmx_dF_cpp(vec({tau_Dmax}), fd32_const1, FDdz_mx, d).col(0);

  vec aug_3 = aug_31 - aug_32;
  vec aug_3_const1 = aug_31_const1 - aug_32_const1;

  return List::create(
    Named("aug_3") = aug_3,
    Named("aug_3_const1") = aug_3_const1
  );
} //  user system elapsed: 32.646   3.165  36.388 






// [[Rcpp::export]]
List truncAIPW_transMean_EF_cpp(int nn,
                                const vec& time,
                                const vec& Q,
                                const mat& Fuz_mx,
                                const mat& Gvz_mx,
                                const vec& u,
                                const vec& v,
                                const vec& nu_time,
                                const vec& nuu,
                                double trim = 0.0) {
  
  double tau2 = v.max() + 1e-10;
  double tau_Tmax = u.max() + 1.0;
  
  // Compute CDF evaluations
  vec Gtz = CDF_eval_mx_cpp(time, Gvz_mx, v).diag();
  vec Gqz = CDF_eval_mx_cpp(Q, Gvz_mx, v).diag();
  vec Fqz = CDF_eval_mx_cpp(Q, Fuz_mx, u).diag();
  
  vec DDen1 = 1.0 / clamp(Gtz, trim, datum::inf);
  vec DDen2 = Fqz / (clamp(Gqz, trim, datum::inf) % clamp(1.0 - Fqz, trim, datum::inf));
  
  mat Fvz_mx = CDF_eval_mx_cpp(v, Fuz_mx, u);
  mat Gvz_trim_sq = pow(clamp(Gvz_mx, trim, datum::inf), 2);
  mat one_minus_Fvz = clamp(1.0 - Fvz_mx, trim, datum::inf);
  
  // Compute at-risk indicator matrix
  mat Q_mat = repmat(Q, 1, v.n_elem);        // nn x |v|
  mat time_mat = repmat(time, 1, v.n_elem);  // nn x |v|
  rowvec v_row = v.t();                      // 1 x |v|
  mat v_mat = repmat(v_row, nn, 1);          // nn x |v|
  
  mat atrisk_mx = conv_to<mat>::from((Q_mat <= v_mat) % (v_mat < time_mat));
  
  
  mat f_mx = atrisk_mx % Fvz_mx / (Gvz_trim_sq % one_minus_Fvz);
  vec DDen3 = int_fmx_dF_cpp(vec({tau2}), f_mx, Gvz_mx, v).col(0);
  
  vec NNum1 = nu_time / clamp(Gtz, trim, datum::inf);
  
  mat nuu_mx = repmat(nuu.t(), nn, 1);
  vec mqz = int_fmx_dF_cpp(Q, nuu_mx, Fuz_mx, u).diag();
  vec NNum2 = mqz / (clamp(Gqz, trim, datum::inf) % clamp(1.0 - Fqz, trim, datum::inf));
  
  mat mvz_mx = int_fmx_dF_cpp(v, nuu_mx, Fuz_mx, u);
  mat fnu_mx = atrisk_mx % mvz_mx / (Gvz_trim_sq % one_minus_Fvz);
  vec NNum3 = int_fmx_dF_cpp(vec({tau2}), fnu_mx, Gvz_mx, v).col(0);
  
  vec Num_AIPW = NNum1 + NNum2 - NNum3;
  vec Den_AIPW = DDen1 + DDen2 - DDen3;
  
  vec DDen4 = 1.0 / clamp(1.0 - Fqz, trim, datum::inf);
  vec NNum4 = nu_time + mqz / clamp(1.0 - Fqz, trim, datum::inf);
  
  vec Enutz = int_fmx_dF_cpp(vec({tau_Tmax}), nuu_mx, Fuz_mx, u).col(0);
  vec NNum5 = Enutz / clamp(1.0 - Fqz, trim, datum::inf);
  
  return List::create(
    Named("Num_AIPW") = Num_AIPW,
    Named("Den_AIPW") = Den_AIPW,
    Named("Num_IPW.Q") = NNum1,
    Named("Den_IPW.Q") = DDen1,
    Named("Num_Reg.T1") = NNum4,
    Named("Num_Reg.T2") = NNum5,
    Named("Den_Reg") = DDen4
  );
}








/*** R
# R-wrappers
aug_QD <- function(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                   u, v, d, Sdyz_vec = NA, trim = 1e-7) {
  
  if (sum(is.na(Sdyz_vec)) > 0) {
    X_res <- X - Q
    Sdyz_vec <- diag(1 - CDF_eval_mx_cpp(X_res, 1 - Sdz_mx, d))
  }
  
  aug_result_1 <- aug_QD_1_cpp1(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                               u = u, v = v, d = d, trim = trim)
  
  aug_result_2 <- aug_QD_2_cpp1(nn, nuu, Q, Delta, Fuz_mx, Gvz_mx,
                               u = u, v = v, Sdyz_vec = Sdyz_vec, trim = trim)
  
  aug_result_3 <- aug_QD_3_cpp1(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                               u = u, v = v, d = d, trim = trim)
  
  Aug_QD_nu <- aug_result_1$aug_1 + aug_result_2$aug_2 - aug_result_3$aug_3
  Aug_QD_const1 <- aug_result_1$aug_1_const1 + aug_result_2$aug_2_const1 - aug_result_3$aug_3_const1
  
  return(list(Aug_QD_nu = as.numeric(Aug_QD_nu),
              Aug_QD_const1 = as.numeric(Aug_QD_const1)))
}

aug_QD_1_cpp <- function(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                         u, v, d, trim){

  result = aug_QD_1_cpp1(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                         u, v, d, trim = trim)

  result <- lapply(result, as.vector)

  return(result)
}

aug_QD_2_cpp <- function(nn, nuu, Q, Delta, Fuz_mx, Gvz_mx,
                         u, v, Sdyz_vec, trim = 1e-7) {

  result <- aug_QD_2_cpp1(nn, nuu, Q, Delta, Fuz_mx, Gvz_mx,
                          u, v, Sdyz_vec, trim = trim)

  result <- lapply(result, as.vector)

  return(result)
}

aug_QD_3_cpp <- function(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                         u, v, d, trim = 1e-7) {

  result <- aug_QD_3_cpp1(nn, nuu, X, Q, Delta, Fuz_mx, Gvz_mx, Sdz_mx,
                          u, v, d, trim = trim)

  result <- lapply(result, as.vector)

  return(result)
}

*/









