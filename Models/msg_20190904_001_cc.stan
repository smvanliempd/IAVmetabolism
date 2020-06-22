// based on msg_20190605_003_cc.stan
// Included calculations of 'average distributions' for areas of the subjects;
//  based on the a_gr parameters while sigma_an = 0.
data{
  int<lower=1> N;
  int<lower=1> N_gr;
  int idx_an[N];
  int idx_gr[N];
  real mu_data_l;
  real sd_data_l;
  real Area_lcs[N];
}
parameters{
  real a;
  vector[20] a_an; 
  vector[16] a_gr;
  real<lower=0> sigma;
  real<lower=0> sigma_an;
}
model{
  vector[N]  mu;
  sigma    ~ cauchy(0 , 1) ;
  sigma_an ~ cauchy(0 , 1);
  a        ~ normal( 0 , 1);
  a_an     ~ normal( 0 , 1);
  a_gr     ~ normal( 0 , 1);
  
  for (i in 1:N) mu[i] = a  +  a_an[idx_an[i]] * sigma_an + a_gr[idx_gr[i]];
  Area_lcs ~ normal(mu, sigma);
}
generated quantities{
  real del_c_3i_0i;
  real del_c_5i_0i;
  real del_c_8i_0i;
  real del_c_18i_0i;
  real del_c_30i_0i;
  real del_c_3i_3m;
  real del_c_5i_5m;
  real del_c_8i_8m;
  real del_c_18i_18m;
  real del_c_30i_30m;
  real del_d_3i_0i;
  real del_d_5i_0i;
  real del_d_3i_3m;
  real del_d_5i_5m;
  real del_d3m_c3m;
  real del_d5m_c5m;
  real del_d0i_c0i;
  real del_d3i_c3i;
  real del_d5i_c5i;
  // real del_del_cd_3i0i;
  vector[N_gr] area_gr_l_hat;
  vector[N_gr] area_gr_hat;
  
  // average subject areas
  for (i in 1:N_gr) area_gr_l_hat[i] = sd_data_l * normal_rng( a + a_gr[i], sigma ) + mu_data_l  ;
  
  // relative differences in infection wrt t=0 in C57
  del_c_3i_0i  = log2(exp(area_gr_l_hat[7] - area_gr_l_hat[6]));
  del_c_5i_0i  = log2(exp(area_gr_l_hat[8] - area_gr_l_hat[6]));
  del_c_8i_0i  = log2(exp(area_gr_l_hat[9] - area_gr_l_hat[6]));
  del_c_18i_0i = log2(exp(area_gr_l_hat[10] - area_gr_l_hat[6]));
  del_c_30i_0i = log2(exp(area_gr_l_hat[11] - area_gr_l_hat[6]));

  // relative differences wrt mock in C57
  del_c_3i_3m   = log2(exp(area_gr_l_hat[7]- area_gr_l_hat[1]));
  del_c_5i_5m   = log2(exp(area_gr_l_hat[8]- area_gr_l_hat[2]));
  del_c_8i_8m   = log2(exp(area_gr_l_hat[9]- area_gr_l_hat[3]));
  del_c_18i_18m = log2(exp(area_gr_l_hat[10] - area_gr_l_hat[4]));
  del_c_30i_30m = log2(exp(area_gr_l_hat[11] - area_gr_l_hat[5]));

  // relative differences wrt t=0 in DBA
  del_d_3i_0i = log2(exp(area_gr_l_hat[15] - area_gr_l_hat[14]));
  del_d_5i_0i = log2(exp(area_gr_l_hat[16] - area_gr_l_hat[14]));
 
 // relative differences wrt mock in DBA
  del_d_3i_3m = log2(exp(area_gr_l_hat[15] - area_gr_l_hat[12]));
  del_d_5i_5m = log2(exp(area_gr_l_hat[16] - area_gr_l_hat[13]));

  // relative inter-strain differences between mock wrt to 
  del_d3m_c3m = log2(exp(area_gr_l_hat[12] - area_gr_l_hat[1]));
  del_d5m_c5m = log2(exp(area_gr_l_hat[13] - area_gr_l_hat[2]));

  // relative inter-strain differences between infected 
  del_d0i_c0i = log2(exp(area_gr_l_hat[14] - area_gr_l_hat[6]));
  del_d3i_c3i = log2(exp(area_gr_l_hat[15] - area_gr_l_hat[7]));
  del_d5i_c5i = log2(exp(area_gr_l_hat[16] - area_gr_l_hat[8]));
  
  // delta-delta c57-dba t0-t3 post-infection
  // del_del_cd_3i0i = log2(exp2(del_d_3i_0i) / exp2(del_c_3i_0i));
  
  // calculate true areas
  for (i in 1:N_gr) area_gr_hat[i] = exp(area_gr_l_hat[i])  ;
  
 //  // relative differences in infection wrt t=0 in C57
 //  del_c_3i_0i  = log2(exp(sd_data_l * (a_gr[7]  - a_gr[6] )) );
 //  del_c_5i_0i  = log2(exp(sd_data_l * (a_gr[8]  - a_gr[6] )) );
 //  del_c_8i_0i  = log2(exp(sd_data_l * (a_gr[9]  - a_gr[6] )) );
 //  del_c_18i_0i = log2(exp(sd_data_l * (a_gr[10] - a_gr[6] )) );
 //  del_c_30i_0i = log2(exp(sd_data_l * (a_gr[11] - a_gr[6] )) );
 // 
 //  // relative differences wrt mock in C57
 //  del_c_3i_3m   = log2(exp(sd_data_l * (a_gr[7]  - a_gr[1] )) );
 //  del_c_5i_5m   = log2(exp(sd_data_l * (a_gr[8]  - a_gr[2] )) );
 //  del_c_8i_8m   = log2(exp(sd_data_l * (a_gr[9]  - a_gr[3] )) );
 //  del_c_18i_18m = log2(exp(sd_data_l * (a_gr[10] - a_gr[4] )) );
 //  del_c_30i_30m = log2(exp(sd_data_l * (a_gr[11] - a_gr[5] )) );
 // 
 //  // relative differences wrt t=0 in DBA
 //  del_d_3i_0i = log2(exp(sd_data_l * (a_gr[15] - a_gr[14] )) );
 //  del_d_5i_0i = log2(exp(sd_data_l * (a_gr[16] - a_gr[14] )) );
 // 
 // // relative differences wrt mock in DBA
 //  del_d_3i_3m = log2(exp(sd_data_l * (a_gr[15] - a_gr[12] )) );
 //  del_d_5i_5m = log2(exp(sd_data_l * (a_gr[16] - a_gr[13] )) );
 // 
 //  // relative inter-strain differences between mock wrt to 
 //  del_d3m_c3m = log2(exp(sd_data_l * (a_gr[12] - a_gr[1] )) );
 //  del_d5m_c5m = log2(exp(sd_data_l * (a_gr[13] - a_gr[2] )) );
 // 
 //  // relative inter-strain differences between infected 
 //  del_d0i_c0i = log2(exp(sd_data_l * (a_gr[14] - a_gr[6] )) );
 //  del_d3i_c3i = log2(exp(sd_data_l * (a_gr[15] - a_gr[7] )) );
 //  del_d5i_c5i = log2(exp(sd_data_l * (a_gr[16] - a_gr[8] )) );
 //  
  
  
}


