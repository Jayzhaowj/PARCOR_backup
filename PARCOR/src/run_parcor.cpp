// PARCOR
// run_parcor is main function of performing PARCOR model
// misc about transformation from PARCOR to AR.
// gen_AR_sample is function of generating sample of AR coefficients.
// PARCOR_to_AR_fun is function of transforming from PARCOR to AR.


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "shared/PARCOR_to_AR.hpp"
#include "shared/mvnorm.h"
#include "shared/gen_F1t.hpp"
//#include "shared/filter.hpp"
#include "shared/DIC.hpp"
//#include <future>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List filter(arma::mat F1_fwd, 
                  arma::mat F1_bwd, 
                  arma::mat G, 
                  arma::mat mk_0,
                  arma::mat Ck_0,
                  double n_0, 
                  arma::mat S_0,
                  int m,
                  arma::rowvec delta,
                  int type_num,
                  int P,
                  int n_t,
                  int n_I,
                  int n_I2
){
  double ll = 0.0;
  int sign = 1;
  arma::dmat at(n_I2, n_t, arma::fill::zeros);
  arma::dmat mt(n_I2, n_t, arma::fill::zeros);
  arma::dmat Rt(n_I2, n_I2, arma::fill::zeros);
  arma::dmat F1t(n_I, n_I2, arma::fill::zeros);
  arma::dcube Ct(n_I2, n_I2, n_t);
  arma::dcube Qt(n_I, n_I, n_t);
  arma::mat St_sqp(n_I, n_I, arma::fill::zeros);
  arma::dcube St(n_I, n_I, n_t);
  arma::dmat S_comp(n_I, n_I, arma::fill::zeros);
  arma::dmat ft(n_I, n_t, arma::fill::zeros);
  arma::dmat At(n_I2, n_I, arma::fill::zeros);
  arma::dmat F1(n_I, n_t, arma::fill::zeros);
  arma::dmat yt(n_I, n_t, arma::fill::zeros);
  int ubound = 0;
  int lbound = 0;
  if(type_num == 1){
    ubound = n_t;
    lbound = P;
    mt.col(P - 1) = mk_0;
    Ct.slice(P - 1) = Ck_0; 
    St.slice(P - 1) = S_0;
    F1 = F1_bwd;
    yt = F1_fwd;
  }else{
    ubound = n_t - P;
    lbound = 0;
    F1 = F1_fwd;
    yt = F1_bwd;
    sign = -1;
  }
  arma::dmat F1_new(n_I, n_t, arma::fill::zeros);
  arma::dmat delta_m = arma::diagmat(arma::pow(delta, -0.5));
  for(int i = lbound; i < ubound; i++){
    F1t = arma::trans(gen_Ft(F1.col(i - sign * m)));
    if(i == 0){
      at.col(i) = G * mk_0;
      Rt = delta_m * G * Ck_0 * arma::trans(G) * delta_m;
      Qt.slice(i) = F1t * Rt * arma::trans(F1t) + S_0;
      St_sqp = arma::sqrtmat_sympd(S_0);
    }else{
      at.col(i) = G * mt.col(i-1);
      Rt = delta_m * G * Ct.slice(i-1) * arma::trans(G) * delta_m;
      Qt.slice(i) = F1t * Rt * arma::trans(F1t) + St.slice(i-1);
      St_sqp = arma::sqrtmat_sympd(St.slice(i-1));
    }
    Rt = 0.5*Rt + 0.5*arma::trans(Rt);
    Qt.slice(i) = 0.5*Qt.slice(i) + 0.5*arma::trans(Qt.slice(i));
    ft.col(i) = F1t * at.col(i);
    arma::dmat Qt_inv = arma::inv_sympd(Qt.slice(i));
    arma::dmat Qt_inv_sq = arma::sqrtmat_sympd(Qt_inv);
    At = Rt * arma::trans(F1t) * Qt_inv;
    arma::colvec et = yt.col(i) - ft.col(i);
    S_comp += St_sqp * Qt_inv_sq * et * arma::trans(et) * Qt_inv_sq * St_sqp;
    St.slice(i) = (n_0*S_0 + S_comp)/(n_0 + i + 1);
    St.slice(i) = 0.5*St.slice(i) + 0.5*arma::trans(St.slice(i));
    Ct.slice(i) = Rt - At * Qt.slice(i) * arma::trans(At);
    mt.col(i) = at.col(i) + At * et;
    if((i >= P) & (i < n_t - P) ){
      arma::vec tmp_ll = dmvnorm(arma::trans(yt.col(i)), ft.col(i), Qt.slice(i), true);
      ll += arma::sum(tmp_ll);
    }
  }
    return Rcpp::List::create(Rcpp::Named("ll") = ll,
                              Rcpp::Named("mt") = mt,
                              Rcpp::Named("at") = at,
                              Rcpp::Named("Ct") = Ct,
                              Rcpp::Named("St") = St,
                              Rcpp::Named("Qt") = Qt,
                              Rcpp::Named("F1") = F1,
                              Rcpp::Named("yt") = yt,
                              Rcpp::Named("ft") = ft);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List filter_smooth(arma::mat F1_fwd, 
                         arma::mat F1_bwd, 
                         arma::mat G,
                         arma::mat mk_0,
                         arma::mat Ck_0,
                         double n_0,
                         arma::mat S_0,
                         int m,
                         arma::mat delta,
                         int type_num,
                         int P,
                         bool DIC,
                         int sample_size,
                         int chains
){
  // initializing
    //std::cout << "a1";
    int n_t = F1_fwd.n_cols;  // the number of time points
    int n_I = F1_fwd.n_rows;  // the number of time series
    int n_I2 = G.n_rows;  // the dimension of each stage
    int delta_n = delta.n_rows; // the number of choice of discount factors
    int lbound = 0;
    int ubound = 0;
    int sign = 1;
    arma::dmat mnt(n_I2, n_t, arma::fill::zeros);
    arma::cube Cnt(n_I2, n_I2, n_t);
    arma::dmat F1(n_I, n_t, arma::fill::zeros);
    arma::dmat yt(n_I, n_t, arma::fill::zeros);
    arma::dmat F1_new(n_I, n_t, arma::fill::zeros);
    arma::rowvec ll(delta_n);
    double ll_DIC = 0.0;
    double es_DIC = 0.0;

    if(type_num == 1){
      F1 = F1_bwd;
      yt = F1_fwd;
    }else{
      F1 = F1_fwd;
      yt = F1_bwd;
    }
    // do the filtering updating function
    Rcpp::List filter_opt = filter(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, m, delta.row(0), type_num, P, n_t, n_I, n_I2);
    double ll_max = filter_opt["ll"];
    arma::rowvec delta_min = delta.row(0);
    ll(0) = ll_max;
    for(int j = 1; j < delta_n; j++){
        Rcpp::List filter_new = filter(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, m, delta.row(j), type_num, P, n_t, n_I, n_I2);
        double ll_new = filter_new["ll"];
        ll(j) = ll_new;
        if(ll_max < ll_new){
            filter_opt = filter_new;
            ll_max = ll_new;
            delta_min = delta.row(j);
        }
      // Rprintf("completation: %i / %i \r", j+1, delta_n);
    }
  
    arma::mat mt = filter_opt["mt"];
    arma::mat at = filter_opt["at"];
    arma::dcube Ct = filter_opt["Ct"];
    arma::mat Rt(n_I2, n_I2, arma::fill::zeros);
    arma::cube St = filter_opt["St"];
  
    if(type_num == 1){
        lbound = P;
        ubound = n_t;
    }else{
        lbound = 0;
        ubound = n_t - P;
        sign = -1;
    }
    // This following part is smoothing.
    // initializing.
    mnt.col(ubound - 1) = mt.col(ubound - 1);
    Cnt.slice(ubound - 1) = Ct.slice(ubound - 1);
    arma::mat delta_m = arma::diagmat(arma::pow(delta_min, -0.5));
    for(int i = (ubound - 2); i > (lbound - 1); i--){
        Rt = delta_m * G * Ct.slice(i) * arma::trans(G) * delta_m;
        Rt = 0.5*Rt + 0.5*arma::trans(Rt);
        arma::mat Rtp1_inv = arma::inv_sympd(Rt);
        arma::mat Bt = Ct.slice(i) * arma::trans(G) * Rtp1_inv;
        mnt.col(i) = mt.col(i) + Bt * (mnt.col(i+1) - at.col(i+1));
        Cnt.slice(i) = Ct.slice(i) - Bt * (Rt - Cnt.slice(i+1))*arma::trans(Bt);
        Cnt.slice(i) = 0.5*Cnt.slice(i) + 0.5*arma::trans(Cnt.slice(i));
    }
    for(int i = lbound; i < ubound; i++){
        arma::mat F1t = arma::trans(gen_Ft(F1.col(i - sign * m)));
        F1_new.col(i) = yt.col(i) - F1t * mnt.col(i); //update the next stage prediction error.
    }
    

    if(DIC){
        Rcpp::List tmp_DIC = compute_DIC(filter_opt, sample_size, m-1, P, chains);
        ll_DIC = tmp_DIC["ll"];
        es_DIC = tmp_DIC["ES_mean"];
    }
    
    
    return Rcpp::List::create(Rcpp::Named("F1new") = F1_new,
                              Rcpp::Named("mnt") = mnt,
                              Rcpp::Named("Cnt") = Cnt,
                              Rcpp::Named("ll_max") = ll_max,
                              Rcpp::Named("delta_min") = delta_min,
                              //Rcpp::Named("filter_opt") = filter_opt,
                              Rcpp::Named("ll") = ll,
                              Rcpp::Named("St") = St,
                              Rcpp::Named("ES_mean") = es_DIC,
                              Rcpp::Named("ll_DIC") = ll_DIC);
    
    
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List  run_parcor(arma::mat F1, 
                       arma::mat G, 
                       arma::mat mk_0,
                       arma::mat Ck_0,
                       double n_0, 
                       arma::mat S_0,
                       arma::cube delta,
                       int P,
                       int sample_size,
                       int chains,
                       bool DIC,
                       bool backward
){
    int n_t = F1.n_cols;
    int n_I = F1.n_rows;
    int n_I2 = G.n_rows;
    int n_delta = delta.slice(0).n_rows;
    arma::mat F1_fwd(n_I, n_t, arma::fill::zeros);
    arma::mat F1_bwd(n_I, n_t, arma::fill::zeros);
    Rcpp::List phi_fwd(P);
    Rcpp::List phi_bwd(P);
    arma::cube akm_prev(n_I, n_t, P);
    arma::cube dkm_prev(n_I, n_t, P);
    arma::vec ll_fwd(P);
    arma::vec ll_bwd(P);
    arma::mat ll(P, n_delta, arma::fill::zeros);
    arma::mat delta_fwd(P, n_I2, arma::fill::zeros);
    arma::mat delta_bwd(P, n_I2, arma::fill::zeros);
    Rcpp::List St_fwd(P);
    Rcpp::List St_bwd(P);
    Rcpp::List ar_coef_prev;
    Rcpp::List ar_coef(P);
    Rcpp::List Cnt_fwd(P);
    Rcpp::List Cnt_bwd(P);
    //Rcpp::List filter_fwd(P);
    //Rcpp::List filter_bwd(P);
    arma::vec val_DIC(P, arma::fill::zeros);
    arma::vec ll_DIC(P, arma::fill::zeros);
    arma::vec es_DIC(P, arma::fill::zeros); //effective sample size
    
  //initialize
    F1_fwd = F1;
    F1_bwd = F1;
    Rcpp::List tmp_fwd;
    Rcpp::List tmp_bwd;
    for(int i = 0; i < P; i++){
        //Rprintf("a \n");
        //std::future<Rcpp::List> t1 = std::async(filter_smooth, F1_bwd, G, F1_fwd, mk_0, Ck_0, n_0, S_0, i+1, delta.slice(i), 1, P, DIC, sample_size, chains);
        //Rprintf("b \n");
        //std::future<Rcpp::List> t2 = std::async(filter_smooth, F1_fwd, G, F1_bwd, mk_0, Ck_0, n_0, S_0, i+1, delta.slice(i), -1, P, false, sample_size, chains);
        
        //Rcpp::List tmp_fwd = t1.get();
        //Rcpp::List tmp_bwd = t2.get();
        tmp_fwd = filter_smooth(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, i+1, delta.slice(i), 1, P, DIC,
                                           sample_size, chains);
        Rprintf("Forward computation is completed!\n");
        delta_fwd.row(i) = Rcpp::as<arma::mat>(tmp_fwd["delta_min"]);
        if(backward){
            tmp_bwd = filter_smooth(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, i+1, delta_fwd.row(i), -1, P, false,
                                    sample_size, chains);
        }else{
            tmp_bwd = filter_smooth(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, i+1, delta.slice(i), -1, P, false,
                                    sample_size, chains);
        }
        Rprintf("Backward computation is completed!\n");
        delta_bwd.row(i) = Rcpp::as<arma::mat>(tmp_bwd["delta_min"]);
        ll_DIC(i) = tmp_fwd["ll_DIC"];
        es_DIC(i) = tmp_fwd["ES_mean"];
        
        Cnt_fwd(i) = tmp_fwd["Cnt"];
        Cnt_bwd(i) = tmp_bwd["Cnt"];
        //filter_fwd(i) = tmp_fwd["filter_opt"];
        //filter_bwd(i) = tmp_bwd["filter_opt"];
        St_fwd(i) = tmp_fwd["St"];
        St_bwd(i) = tmp_bwd["St"];
        F1_fwd = Rcpp::as<arma::mat>(tmp_fwd["F1new"]);
        F1_bwd = Rcpp::as<arma::mat>(tmp_bwd["F1new"]);
        phi_fwd(i) = Rcpp::as<arma::mat>(tmp_fwd["mnt"]);
        phi_bwd(i) = Rcpp::as<arma::mat>(tmp_bwd["mnt"]);
        ll_fwd(i) = tmp_fwd["ll_max"];
        ll_bwd(i) = tmp_bwd["ll_max"];
        ll.row(i) = Rcpp::as<arma::rowvec>(tmp_fwd["ll"]);
        
        
        if(i == 0){
            ar_coef(i) = PARCOR_to_AR(Rcpp::as<arma::mat>(phi_fwd(i)), Rcpp::as<arma::mat>(phi_bwd(i)), akm_prev, dkm_prev, n_I, i, P);
        }else{
            ar_coef_prev = ar_coef(i - 1);
            akm_prev = Rcpp::as<arma::cube>(ar_coef_prev["forward"]);
            dkm_prev = Rcpp::as<arma::cube>(ar_coef_prev["backward"]);
            ar_coef(i) = PARCOR_to_AR(Rcpp::as<arma::mat>(phi_fwd(i)), Rcpp::as<arma::mat>(phi_bwd(i)), akm_prev, dkm_prev, n_I, i, P);
        }
        Rprintf("iterations: %i / %i \n", i+1,  P);
        Rprintf("Forward likelihood at %i stage: %f; Backward likelihood at %i stage: %f. \n", i+1, ll_fwd(i), i+1, ll_bwd(i));
    }
    val_DIC = 2*(arma::cumsum(es_DIC) - ll_DIC);
    return Rcpp::List::create(Rcpp::Named("phi_fwd") = phi_fwd,
                              Rcpp::Named("phi_bwd") = phi_bwd,
                              Rcpp::Named("St_fwd") = St_fwd,
                              Rcpp::Named("St_bwd") = St_bwd,
                              Rcpp::Named("ll_fwd") = ll_fwd,
                              Rcpp::Named("ll_bwd") = ll_bwd,
                              Rcpp::Named("delta_fwd") = delta_fwd,
                              Rcpp::Named("delta_bwd") = delta_bwd,
                              Rcpp::Named("ll") = ll,
                              Rcpp::Named("ar_coef") = ar_coef,
                              Rcpp::Named("Cnt_fwd") = Cnt_fwd,
                              Rcpp::Named("Cnt_bwd") = Cnt_bwd,
                              Rcpp::Named("DIC") = val_DIC,
                              Rcpp::Named("ES_mean") = es_DIC,
                              Rcpp::Named("ll_DIC") = ll_DIC);
                              //Rcpp::Named("Cnt_fwd") = Cnt_fwd,
                              //Rcpp::Named("Cnt_bwd") = Cnt_bwd,
                              //Rcpp::Named("filter_fwd") = filter_fwd,
                              //Rcpp::Named("filter_bwd") = filter_bwd,);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List  gen_AR_sample(arma::cube phi_fwd,
                          arma::cube phi_bwd,
                          Rcpp::List Cnt_fwd,
                          Rcpp::List Cnt_bwd,
                          int n_I,
                          int P_opt,
                          int P_max,
                          int h
                          ){
    int n_t = phi_fwd.n_cols;
    int n_I2 = n_I * n_I;
    arma::cube akm_prev(n_I, n_t, P_opt);
    arma::cube dkm_prev(n_I, n_t, P_opt);
    Rcpp::List ar_coef_prev;
    Rcpp::List ar_coef(P_opt);
    //initialize
    for(int i = 0; i < P_opt; i++){
      arma::mat phi_fwd_sample(n_I2, n_t, arma::fill::zeros);
      arma::mat phi_bwd_sample(n_I2, n_t, arma::fill::zeros);
      arma::mat mnt_fwd = phi_fwd.slice(i);
      arma::mat mnt_bwd = phi_bwd.slice(i);
      arma::cube Cnt_fwd_cur = Cnt_fwd(i);
      arma::cube Cnt_bwd_cur = Cnt_bwd(i);
      for(int j = n_t-1; j > (P_max-1); j--){

          phi_fwd_sample.col(j) = arma::trans(rmvnorm(1, mnt_fwd.col(j), Cnt_fwd_cur.slice(j)));
          if(j > (n_t-h-P_max-1)){
            phi_bwd_sample.col(j) = phi_fwd_sample.col(j);
          }else{
            phi_bwd_sample.col(j) = arma::trans(rmvnorm(1, mnt_bwd.col(j), Cnt_bwd_cur.slice(j)));
          }
      }
        if(i == 0){
            ar_coef(i) = PARCOR_to_AR(phi_fwd_sample, phi_bwd_sample, akm_prev, dkm_prev, n_I, i, P_max);
        }else{
            ar_coef_prev = ar_coef(i - 1);
            akm_prev = Rcpp::as<arma::cube>(ar_coef_prev["forward"]);
            dkm_prev = Rcpp::as<arma::cube>(ar_coef_prev["backward"]);
            ar_coef(i) = PARCOR_to_AR(phi_fwd_sample, phi_bwd_sample, akm_prev, dkm_prev, n_I, i, P_max);
        }
    }
    return ar_coef;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List PAR_to_AR_fun(arma::cube phi_fwd, arma::cube phi_bwd, int n_I){
    int n_t = phi_fwd.n_cols;
    int P = phi_fwd.n_slices;
    arma::cube akm_prev(n_I, n_t, P);
    arma::cube dkm_prev(n_I, n_t, P);
    Rcpp::List ar_coef_prev;
    Rcpp::List ar_coef(P);
    for(int i = 0; i < P; i++){
        if(i == 0){
            ar_coef(i) = PARCOR_to_AR(phi_fwd.slice(i), phi_bwd.slice(i), akm_prev, dkm_prev, n_I, i, P);
        }else{
            ar_coef_prev = ar_coef(i - 1);
            akm_prev = Rcpp::as<arma::cube>(ar_coef_prev["forward"]);
            dkm_prev = Rcpp::as<arma::cube>(ar_coef_prev["backward"]);
            ar_coef(i) = PARCOR_to_AR(phi_fwd.slice(i), phi_bwd.slice(i), akm_prev, dkm_prev, n_I, i, P);
        }
    }
    return ar_coef;
}


