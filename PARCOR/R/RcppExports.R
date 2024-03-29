# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cp_sd <- function(phi, SIGMA, w) {
    .Call(`_PARCOR_cp_sd`, phi, SIGMA, w)
}

get_sd <- function(sd, ts1, ts2, type) {
    .Call(`_PARCOR_get_sd`, sd, ts1, ts2, type)
}

compute_DIC_TVVAR <- function(temp_filter, sample_size, St, P_max) {
    .Call(`_PARCOR_compute_DIC_TVVAR`, temp_filter, sample_size, St, P_max)
}

compute_spec <- function(phi, SIGMA, w, P_max, ch1, ch2, time_depend = TRUE) {
    .Call(`_PARCOR_compute_spec`, phi, SIGMA, w, P_max, ch1, ch2, time_depend)
}

filter_smooth_TVVAR <- function(F1, G, mk_0, Ck_0, n_0, S_0, m, delta, pp) {
    .Call(`_PARCOR_filter_smooth_TVVAR`, F1, G, mk_0, Ck_0, n_0, S_0, m, delta, pp)
}

filter <- function(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, m, delta, type_num, P, n_t, n_I, n_I2) {
    .Call(`_PARCOR_filter`, F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, m, delta, type_num, P, n_t, n_I, n_I2)
}

filter_smooth <- function(F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, m, delta, type_num, P, DIC, sample_size, chains, uncertainty) {
    .Call(`_PARCOR_filter_smooth`, F1_fwd, F1_bwd, G, mk_0, Ck_0, n_0, S_0, m, delta, type_num, P, DIC, sample_size, chains, uncertainty)
}

run_parcor <- function(F1, G, mk_0, Ck_0, n_0, S_0, delta, P, sample_size, chains, DIC, backward, uncertainty) {
    .Call(`_PARCOR_run_parcor`, F1, G, mk_0, Ck_0, n_0, S_0, delta, P, sample_size, chains, DIC, backward, uncertainty)
}

gen_AR_sample <- function(phi_fwd, phi_bwd, Cnt_fwd, Cnt_bwd, n_I, P_opt, P_max, h) {
    .Call(`_PARCOR_gen_AR_sample`, phi_fwd, phi_bwd, Cnt_fwd, Cnt_bwd, n_I, P_opt, P_max, h)
}

PAR_to_AR_fun <- function(phi_fwd, phi_bwd, n_I) {
    .Call(`_PARCOR_PAR_to_AR_fun`, phi_fwd, phi_bwd, n_I)
}

