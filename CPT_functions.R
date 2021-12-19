# CPT estimation with R and stan
# Author: Alina Ferecatu
# Date: 12/15/2021

## Basic functions ----
logsumexp <- function(x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function(x) {
  exp(x - logsumexp(x))
}

inv_logit <- function(logit) {
  exp(logit) / (1 + exp(logit))
}

### PWF Gonzalez Wu 1999 ----
pwf_gonzalezwu<-function(p, delta, gamma){
  wp=(delta*(p^gamma))/(delta*(p^gamma)+((1-p)^gamma))
  return(wp)
}

### stan code CPT aggregate model for ladders ----
CPT_agg_ladders<-"
functions {
//gonzalez wu 1999;
real pwf_gonzalezwu(real p, real delta, real gamma){
real wp;
wp = (delta*(p^gamma))/(delta*(p^gamma)+((1-p)^gamma));
return wp;
}

}

data {
int Nobs; // number of obs;
int Nind; // number of subject;
int id[Nobs];
int<lower=1, upper=2> decision[Nobs]; //1-risky, 2-safe;

real<lower=0, upper=1> opt1prob1[Nobs];
real<lower=0, upper=1> opt1prob2[Nobs];
real<lower=0, upper=1> opt2prob1[Nobs];

real opt1val1[Nobs];
real opt1val2[Nobs];
real opt2val1[Nobs];

int<lower=1> nvar;
}

parameters{
//group-level parameters only
vector[nvar] beta; //
}

transformed parameters {
//subject-level parameters
real<lower=0, upper=1> gamma;
real<lower=0> delta;
real<lower=0, upper=2> alpha;
real<lower=0> phi; // sensitivity parameter

// transform parameter to stay within appropriate ranges;
gamma = Phi_approx(beta[1]);
delta = exp(beta[2]);
alpha  = Phi_approx(beta[3])*2;
phi    = exp(beta[4]);
}

model {
// draw untransformed parameters;
to_vector(beta) ~ normal(0, 1);

//subject loop and trial loop;
for (n in 1:Nobs) {
vector[3] w_prob;
vector[2] U_opt;

//probability weight function - change funtion here;
  w_prob[1] = pwf_gonzalezwu(opt1prob1[n], delta, gamma);
  // CPT - pi(p+q)-pi(p);
  w_prob[2] = pwf_gonzalezwu((opt1prob1[n]+opt1prob2[n]), delta, gamma)-w_prob[1];
  w_prob[3] = pwf_gonzalezwu(opt2prob1[n], delta, gamma); 
 
//compute subjective value of option 1;
U_opt[1] = w_prob[1]*(opt1val1[n]^alpha ) + w_prob[2]*(opt1val2[n]^alpha );

//compute subjective value of option 2;
U_opt[2]  = w_prob[3]*(opt2val1[n]^alpha );

// compute action probabilities
decision[n] ~ categorical_logit(U_opt*phi );
}
}

generated quantities {

real log_lik[Nobs];
// For posterior predictive check
real y_pred[Nobs];
real pred_prob[Nobs];

{ // local section, this saves time and space

for (n in 1:Nobs) {

vector[3] w_prob;
vector[2] U_opt;

//probability weight function - change funtion here;
w_prob[1] = pwf_gonzalezwu(opt1prob1[n], delta, gamma);
//CPT - pi(p+q)-pi(p);
w_prob[2] = pwf_gonzalezwu((opt1prob1[n]+opt1prob2[n]), delta, gamma)-w_prob[1];
w_prob[3] = pwf_gonzalezwu(opt2prob1[n], delta, gamma); 
 
//compute subjective value of option 1;
U_opt[1] = w_prob[1]*(opt1val1[n]^alpha ) + w_prob[2]*(opt1val2[n]^alpha );

//compute subjective value of option 2;
U_opt[2]  = w_prob[3]*(opt2val1[n]^alpha );

// compute action probabilities
log_lik[n]=categorical_logit_lpmf(decision[n]|U_opt*phi );
y_pred[n]  = categorical_rng(softmax(U_opt*phi));
pred_prob[n]= softmax(U_opt*phi )[1];
}
}
}
"


### stan code CPT  individual model (with correlated priors) for ladders ----
CPT_indcorr_ladders <- "
functions {
//gonzalez wu 1999;
real pwf_gonzalezwu(real p, real delta, real gamma){
real wp;
wp = (delta*(p^gamma))/(delta*(p^gamma)+((1-p)^gamma));
return wp;
}

}

data {
int Nobs; // number of obs;
int Nind; // number of subject;
int id[Nobs];
int<lower=1, upper=2> decision[Nobs]; //1-risky, 2 safe;
real<lower=0, upper=1> opt1prob1[Nobs];
real<lower=0, upper=1> opt1prob2[Nobs];
real<lower=0, upper=1> opt2prob1[Nobs];

real opt1val1[Nobs];
real opt1val2[Nobs];
real opt2val1[Nobs];
int<lower=1> nvar;

int<lower=1> ndem;
matrix[Nind, ndem] u;

}

parameters{
//group-level parameters
matrix[nvar, Nind] nu; // nvar*Nind parameter matrix
matrix[ndem, nvar] tau; // 
cholesky_factor_corr[nvar] L_Omega; // Vbeta - prior correlation - off diagonal elements
vector<lower=0>[nvar] sigma; //standard deviation - diagonal elements
}

transformed parameters {
//subject-level parameters
real<lower=0, upper=1> gamma[Nind];
real<lower=0> delta[Nind];
real<lower=0, upper=2> alpha[Nind];
real<lower=0> phi[Nind];

// noncentered reparametrisation
matrix[Nind, nvar] beta;
matrix[Nind, nvar] Vbeta_reparametrized;

Vbeta_reparametrized = (diag_pre_multiply(sigma, L_Omega)*nu)';
beta=u*tau+Vbeta_reparametrized;

for (i in 1:Nind) {
gamma[i] = Phi_approx(beta[i, 1]);
delta[i] = exp(beta[i, 2]);
alpha[i] = Phi_approx(beta[i, 3])*2;
phi[i]   = exp(beta[i, 4]);
}

}

model {
//prior : hyperparameters
L_Omega~lkj_corr_cholesky(2);
sigma~exponential(1);
to_vector(nu)~normal(0, 1);
to_vector(tau) ~ normal(0, 1);

//subject loop and trial loop;
for (n in 1:Nobs) {
vector[3] w_prob;
vector[2] U_opt;

//probability weight function - change funtion here;
  w_prob[1] = pwf_gonzalezwu(opt1prob1[n], delta[id[n]], gamma[id[n]]);
  //pi(p+q)-pi(p);
  w_prob[2] = pwf_gonzalezwu((opt1prob1[n]+opt1prob2[n]), delta[id[n]], gamma[id[n]])-w_prob[1];
  w_prob[3] = pwf_gonzalezwu(opt2prob1[n], delta[id[n]], gamma[id[n]]);
  
//compute subjective value of option 1;
U_opt[1] = w_prob[1]*(opt1val1[n]^alpha[id[n]]) + w_prob[2]*(opt1val2[n]^alpha[id[n]]);

//compute subjective value of option 2;
U_opt[2]  = w_prob[3]*(opt2val1[n]^alpha[id[n]]); 

//print(U_opt);
// compute action probabilities
decision[n] ~ categorical_logit(U_opt*phi[id[n]]);
}
}

generated quantities {
real<lower = 0, upper = 1> mu_gamma;
real<lower = 0> mu_delta;
real<lower = 0, upper = 2> mu_alpha;
real<lower = 0> mu_phi;

real log_lik[Nind];
// For posterior predictive check
real y_pred[Nobs];
real pred_prob[Nobs];
corr_matrix[nvar] Omega;

mu_gamma = Phi_approx( tau[1,1]);
mu_delta = exp( tau[1,2]);
mu_alpha  = Phi_approx(tau[1,3])*2;
mu_phi    = exp(tau[1,4]);

Omega = L_Omega * L_Omega';

{ // local section, this saves time and space
//initialize loglik
for (i in 1:Nind){
log_lik[i] = 0;}

for (n in 1:Nobs) {

vector[3] w_prob;
vector[2] U_opt;

//probability weight function - change funtion here;
  w_prob[1] = pwf_gonzalezwu(opt1prob1[n], delta[id[n]], gamma[id[n]]);
  //CPT - pi(p+q)-pi(p);
  w_prob[2] = pwf_gonzalezwu((opt1prob1[n]+opt1prob2[n]), delta[id[n]], gamma[id[n]])-w_prob[1];
  w_prob[3] = pwf_gonzalezwu(opt2prob1[n], delta[id[n]], gamma[id[n]]);
  
 
//compute subjective value of option 1;
U_opt[1] = w_prob[1]*(opt1val1[n]^alpha[id[n]]) + w_prob[2]*(opt1val2[n]^alpha[id[n]]);

//compute subjective value of option 2;
U_opt[2]  = w_prob[3]*(opt2val1[n]^alpha[id[n]]); 

// compute action probabilities
log_lik[id[n]] = log_lik[id[n]]+categorical_logit_lpmf(decision[n]|U_opt*phi[id[n]]);
y_pred[n]  = categorical_rng(softmax(U_opt*phi[id[n]]));
pred_prob[n]= softmax(U_opt*phi[id[n]])[1];
}
}
}
"


### stan code CPT aggregate model for two-outcome gamble pairs ----
CPT_agg_two_outcome_gambles<-"
functions {
//gonzalez wu 1999;
real pwf_gonzalezwu(real p, real delta, real gamma){
real wp;
wp = (delta*(p^gamma))/(delta*(p^gamma)+((1-p)^gamma));
return wp;
}
}

data {
int Nobs; // number of subjects;
int Nind; // number of obs per subject;
int id[Nobs];
int<lower=1, upper=2> decision[Nobs];
real<lower=0, upper=1> opt1prob[Nobs];
real<lower=0, upper=1> opt2prob[Nobs];
real opt1val1[Nobs];
real opt1val2[Nobs];
real opt2val1[Nobs];
real opt2val2[Nobs];

int<lower=1> nvar;
}

parameters{
//group-level parameters only
vector[nvar] beta; 
}

transformed parameters {

real<lower=0, upper=1> gamma;
real<lower=0> delta;
real<lower=0, upper=2> alpha; 
real<lower=0, upper=5> lambda; 
real<lower=0> phi;

gamma  = Phi_approx(beta[1]);
delta  = exp(beta[2]);
alpha  = Phi_approx(beta[3])*2;
lambda = Phi_approx(beta[4])*5;
phi    = exp(beta[5]);

}

model {

to_vector(beta) ~ normal(0, 1);

//subject loop and trial loop;
for (n in 1:Nobs) {
vector[4] w_prob;
vector[2] U_opt;
vector[4] U_val;

//probability weight function 
if (opt1val1[n]>= 0) {w_prob[1] = pwf_gonzalezwu(opt1prob[n],  delta, gamma);
  if (opt1val2[n]>= 0) w_prob[2] = 1-w_prob[1];
      else  w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta, gamma);
}else{ w_prob[1] = pwf_gonzalezwu(opt1prob[n], delta, gamma);
  if (opt1val2[n]>= 0) w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta, gamma);
    else w_prob[2] = 1-w_prob[1];
}
    
if (opt2val1[n]>= 0) {w_prob[3] = pwf_gonzalezwu(opt2prob[n],  delta, gamma);
  if (opt2val2[n] >= 0) w_prob[4] = 1-w_prob[3];
      else w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta, gamma);
}else{ w_prob[3] = pwf_gonzalezwu(opt2prob[n], delta,  gamma);
  if (opt1val2[n]>= 0) w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta, gamma);
    else w_prob[4] = 1-w_prob[3];
}

//compute subjective value of option 1;
U_val[1] = (opt1val1[n]>= 0 ? w_prob[1]*(opt1val1[n]^alpha) : -w_prob[1]*(fabs(opt1val1[n])^alpha)*lambda );
U_val[2] = (opt1val2[n]>= 0 ? w_prob[2]*(opt1val2[n]^alpha) : -w_prob[2]*(fabs(opt1val2[n])^alpha)*lambda );

//compute subjective value of option 2;
U_val[3] = (opt2val1[n] >= 0 ? w_prob[3]*(opt2val1[n]^alpha) : -w_prob[3]*(fabs(opt2val1[n])^alpha)*lambda ); 
U_val[4] = (opt2val2[n] >= 0 ? w_prob[4]*(opt2val2[n]^alpha) : -w_prob[4]*(fabs(opt2val2[n])^alpha)*lambda );

U_opt[1]=U_val[1]+U_val[2];
U_opt[2]=U_val[3]+U_val[4];

// compute action probabilities
decision[n] ~ categorical_logit(U_opt*phi );
}
}

generated quantities {

real log_lik[Nobs];
// For posterior predictive check
real y_pred[Nobs];
real pred_prob[Nobs];

{ // local section, this saves time and space

for (n in 1:Nobs) {
vector[4] w_prob;
vector[2] U_opt;
vector[4] U_val;

//probability weight function 
if (opt1val1[n]>= 0) {w_prob[1] = pwf_gonzalezwu(opt1prob[n],  delta, gamma);
  if (opt1val2[n]>= 0) w_prob[2] = 1-w_prob[1];
      else  w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta, gamma);
}else{ w_prob[1] = pwf_gonzalezwu(opt1prob[n], delta, gamma);
  if (opt1val2[n]>= 0) w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta, gamma);
    else w_prob[2] = 1-w_prob[1];
}
    
if (opt2val1[n]>= 0) {w_prob[3] = pwf_gonzalezwu(opt2prob[n],  delta, gamma);
  if (opt2val2[n] >= 0) w_prob[4] = 1-w_prob[3];
      else w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta, gamma);
}else{ w_prob[3] = pwf_gonzalezwu(opt2prob[n], delta,  gamma);
  if (opt1val2[n]>= 0) w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta, gamma);
    else w_prob[4] = 1-w_prob[3];
}

//compute subjective value of option 1;
U_val[1] = (opt1val1[n]>= 0 ? w_prob[1]*(opt1val1[n]^alpha) : -w_prob[1]*(fabs(opt1val1[n])^alpha)*lambda );
U_val[2] = (opt1val2[n]>= 0 ? w_prob[2]*(opt1val2[n]^alpha) : -w_prob[2]*(fabs(opt1val2[n])^alpha)*lambda );

//compute subjective value of option 2;
U_val[3] = (opt2val1[n] >= 0 ? w_prob[3]*(opt2val1[n]^alpha) : -w_prob[3]*(fabs(opt2val1[n])^alpha)*lambda ); 
U_val[4] = (opt2val2[n] >= 0 ? w_prob[4]*(opt2val2[n]^alpha) : -w_prob[4]*(fabs(opt2val2[n])^alpha)*lambda );

U_opt[1]=U_val[1]+U_val[2];
U_opt[2]=U_val[3]+U_val[4];

// compute action probabilities
log_lik[n]=categorical_logit_lpmf(decision[n]|U_opt*phi );
y_pred[n]  = categorical_rng(softmax(U_opt*phi));
pred_prob[n]= softmax(U_opt*phi)[1];
}
}
}
"


### stan code CPT individual model (with correlated priors) for two-outcome gamble pairs ----
CPT_indcorr_two_outcome_gambles<-"
functions {
//gonzalez wu 1999;
real pwf_gonzalezwu(real p, real delta, real gamma){
real wp;
wp = (delta*(p^gamma))/(delta*(p^gamma)+((1-p)^gamma));
return wp;
}
}


data {
int Nobs; // number of obs;
int Nind; // number of subject;
int id[Nobs];
int<lower=1, upper=2> decision[Nobs];
real<lower=0, upper=1> opt1prob[Nobs];
real<lower=0, upper=1> opt2prob[Nobs];
real opt1val1[Nobs];
real opt1val2[Nobs];
real opt2val1[Nobs];
real opt2val2[Nobs];

int<lower=1> nvar;
int<lower=1> ndem;
matrix[Nind, ndem] u;
}

parameters{
//group-level parameters
matrix[nvar, Nind] nu; // nvar*Nind parameter matrix
matrix[ndem, nvar] tau; // 
cholesky_factor_corr[nvar] L_Omega; // Vbeta - prior correlation - off-diagonal elements
vector<lower=0>[nvar] sigma; //standard deviation - diagonal elements
}

transformed parameters {
//subject-level parameters
real<lower=0, upper=1> gamma[Nind];
real<lower=0> delta[Nind];
real<lower=0, upper=2> alpha[Nind];
real<lower=0, upper=5> lambda[Nind]; 
real<lower=0> phi[Nind];

// non-centered reparametrisation;
matrix[Nind, nvar] beta;
matrix[Nind, nvar] Vbeta_reparametrized;

Vbeta_reparametrized = (diag_pre_multiply(sigma, L_Omega)*nu)';
beta=u*tau+Vbeta_reparametrized;

for (i in 1:Nind) {
gamma[i] = Phi_approx(beta[i, 1]);
delta[i] = exp(beta[i, 2]);
alpha[i] =  Phi_approx(beta[i, 3])*2;
lambda[i]  = Phi_approx(beta[i, 4])*5;
phi[i]   = exp(beta[i, 5]);
}

}

model {
//prior : hyperparameters
L_Omega~lkj_corr_cholesky(2);
sigma~exponential(1);
to_vector(nu)~normal(0, 1);
to_vector(tau) ~ normal(0, 1);

//subject loop and trial loop;
for (n in 1:Nobs) {

vector[4] w_prob;
vector[2] U_opt;
vector[4] U_val;

//probability weight function 
if (opt1val1[n]>= 0) {w_prob[1] = pwf_gonzalezwu(opt1prob[n],  delta[id[n]], gamma[id[n]]);
  if (opt1val2[n]>= 0) w_prob[2] = 1-w_prob[1];
      else  w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta[id[n]], gamma[id[n]]);
}else{ w_prob[1] = pwf_gonzalezwu(opt1prob[n], delta[id[n]], gamma[id[n]]);
  if (opt1val2[n]>= 0) w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta[id[n]], gamma[id[n]]);
    else w_prob[2] = 1-w_prob[1];
}
    
if (opt2val1[n]>= 0) {w_prob[3] = pwf_gonzalezwu(opt2prob[n],  delta[id[n]], gamma[id[n]]);
  if (opt2val2[n] >= 0) w_prob[4] = 1-w_prob[3];
      else w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta[id[n]], gamma[id[n]]);
}else{ w_prob[3] = pwf_gonzalezwu(opt2prob[n], delta[id[n]],  gamma[id[n]]);
  if (opt1val2[n]>= 0) w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta[id[n]], gamma[id[n]]);
    else w_prob[4] = 1-w_prob[3];
}

//compute subjective value of option 1;
U_val[1] = (opt1val1[n]>= 0 ? w_prob[1]*(opt1val1[n]^alpha[id[n]]) : -w_prob[1]*(fabs(opt1val1[n])^alpha[id[n]] )*lambda[id[n]] );
U_val[2] = (opt1val2[n]>= 0 ? w_prob[2]*(opt1val2[n]^alpha[id[n]]) : -w_prob[2]*(fabs(opt1val2[n])^alpha[id[n]] )*lambda[id[n]] );

//compute subjective value of option 2;
U_val[3] = (opt2val1[n] >= 0 ? w_prob[3]*(opt2val1[n]^alpha[id[n]]) : -w_prob[3]*(fabs(opt2val1[n])^alpha[id[n]] )*lambda[id[n]]); 
U_val[4] = (opt2val2[n] >= 0 ? w_prob[4]*(opt2val2[n]^alpha[id[n]]) : -w_prob[4]*(fabs(opt2val2[n])^alpha[id[n]] )*lambda[id[n]]);

U_opt[1]=U_val[1]+U_val[2];
U_opt[2]=U_val[3]+U_val[4];

// compute action probabilities
decision[n] ~ categorical_logit(U_opt*phi[id[n]] );
}
}

generated quantities {
real<lower = 0, upper = 1> mu_gamma;
real<lower = 0> mu_delta;
real<lower = 0, upper = 2> mu_alpha;
real<lower = 0, upper = 5> mu_lambda;
real<lower = 0> mu_phi;

real log_lik[Nind];
// For posterior predictive check
real y_pred[Nobs];
real pred_prob[Nobs];
corr_matrix[nvar] Omega;

mu_gamma = Phi_approx(tau[1,1]);
mu_delta = exp(tau[1,2]);
mu_alpha  = Phi_approx( tau[1,3])*2;
mu_lambda  = Phi_approx(tau[1,4])*5;
mu_phi    = exp(tau[1,5]);

Omega = L_Omega * L_Omega';

{ // local section, this saves time and space
//initialize loglik
for (i in 1:Nind){
log_lik[i] = 0;}

for (n in 1:Nobs) {

vector[4] w_prob;
vector[2] U_opt;
vector[4] U_val;

//probability weight function 
if (opt1val1[n]>= 0) {w_prob[1] = pwf_gonzalezwu(opt1prob[n],  delta[id[n]], gamma[id[n]]);
  if (opt1val2[n]>= 0) w_prob[2] = 1-w_prob[1];
      else  w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta[id[n]], gamma[id[n]]);
}else{ w_prob[1] = pwf_gonzalezwu(opt1prob[n], delta[id[n]], gamma[id[n]]);
  if (opt1val2[n]>= 0) w_prob[2] = pwf_gonzalezwu((1-opt1prob[n]), delta[id[n]], gamma[id[n]]);
    else w_prob[2] = 1-w_prob[1];
}
    
if (opt2val1[n]>= 0) {w_prob[3] = pwf_gonzalezwu(opt2prob[n],  delta[id[n]], gamma[id[n]]);
  if (opt2val2[n] >= 0) w_prob[4] = 1-w_prob[3];
      else w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta[id[n]], gamma[id[n]]);
}else{ w_prob[3] = pwf_gonzalezwu(opt2prob[n], delta[id[n]],  gamma[id[n]]);
  if (opt1val2[n]>= 0) w_prob[4] = pwf_gonzalezwu((1-opt2prob[n]), delta[id[n]], gamma[id[n]]);
    else w_prob[4] = 1-w_prob[3];
}

//compute subjective value of option 1;
U_val[1] = (opt1val1[n]>= 0 ? w_prob[1]*(opt1val1[n]^alpha[id[n]]) : -w_prob[1]*(fabs(opt1val1[n])^alpha[id[n]] )*lambda[id[n]] );
U_val[2] = (opt1val2[n]>= 0 ? w_prob[2]*(opt1val2[n]^alpha[id[n]]) : -w_prob[2]*(fabs(opt1val2[n])^alpha[id[n]] )*lambda[id[n]] );

//compute subjective value of option 2;
U_val[3] = (opt2val1[n] >= 0 ? w_prob[3]*(opt2val1[n]^alpha[id[n]]) : -w_prob[3]*(fabs(opt2val1[n])^alpha[id[n]] )*lambda[id[n]]); 
U_val[4] = (opt2val2[n] >= 0 ? w_prob[4]*(opt2val2[n]^alpha[id[n]]) : -w_prob[4]*(fabs(opt2val2[n])^alpha[id[n]] )*lambda[id[n]]);

U_opt[1]=U_val[1]+U_val[2];
U_opt[2]=U_val[3]+U_val[4];

// compute action probabilities
log_lik[id[n]] = log_lik[id[n]]+categorical_logit_lpmf(decision[n]|U_opt*phi[id[n]]);
y_pred[n]  = categorical_rng(softmax(U_opt*phi[id[n]]));
pred_prob[n]= softmax(U_opt*phi[id[n]])[1];
}
}
}
" 

in_sample_stats<-function(data, model_stanfit) {
  
  load(paste0(PATH_RESULTS, "/CPT_", model_stanfit, ".RData"))
  
  waic=waic(extract_log_lik(CPT_stanfit))[["estimates"]][3]
  ll= extract_log_lik(CPT_stanfit)
  lpd=waic(ll)[["estimates"]][1] + waic(ll)[["estimates"]][2]
  
  hr_obs=c()
  yrep=as.matrix(CPT_stanfit, pars = c("y_pred"))
  for (j in 1:ncol(yrep)){
    ## hr
    hr_obs[j]=sum(ifelse(yrep[,j]==data$decision[j], 1,0))/nrow(yrep)
  }
  hr=mean(hr_obs)
  
  return(c(round(lpd, 0), round(waic, 0), round(hr, 3) ) )
}



