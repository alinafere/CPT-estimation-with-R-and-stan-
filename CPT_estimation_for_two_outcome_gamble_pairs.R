# CPT estimation with R and stan
# Author: Alina Ferecatu
# Date: 12/15/2021

## Load required packages ----
rm(list=ls())
library(loo)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Set paths and flags ----
## IND =T individul-level model; IND = F aggregate level model
IND=T
PATH_FUNCTION="~/Documents/CPT_estimation_Rstan"
PATH_DATA="~/Documents/CPT_estimation_Rstan/data"
PATH_RESULTS="~/Documents/CPT_estimation_Rstan/results"
PATH_PLOTS= "~/Documents/CPT_estimation_Rstan/plots"

## Source functions ----
source(paste0(PATH_FUNCTION,"/CPT_functions.R") ) 
##* compile stan model (individual or aggregate) ----
if(IND){
  m = stan_model(model_code = CPT_indcorr_two_outcome_gambles)
}else{  m = stan_model(model_code =CPT_agg_two_outcome_gambles)
}  

# Import data ----
rieskamp_gambles=read_csv(paste(PATH_DATA,"/rieskamp08_gambles.csv", sep=""))
rieskamp_choices=read_csv(paste(PATH_DATA,"/rieskamp08_choices.csv", sep=""))
## 0=choice of A (first option); 1=choice of B (second opiton)
rieskamp08_data=rieskamp_choices %>% 
  pivot_longer(-choicepair, names_to="subj_id", values_to="choice") %>% 
  left_join(rieskamp_gambles)%>%
  mutate(decision=(choice+1), subj_id=as.numeric(subj_id)) %>% 
  arrange(subj_id, choicepair)

##* summary stats ----
dim(rieskamp08_data) ## 5400 observations
length(unique(rieskamp08_data$choicepair)) ## 180 unique gambles
length(unique(rieskamp08_data$subj_id)) ## 30 individuals

## choice of options with gambles in the gains domain
table(rieskamp08_data %>% filter(choicepair %in% 1:60) %>% 
        dplyr::select(decision) )
## choice of options with gambles in the loss domain
table(rieskamp08_data %>% filter(choicepair %in% 61:120) %>% 
        dplyr::select(decision) )
## choice of options with mixed gambles
table(rieskamp08_data %>% filter(choicepair %in% 121:180) %>% 
        dplyr::select(decision) )

#* order values to apply CPT probability weighting ----
rieskamp_data=rieskamp08_data %>% 
  mutate(opt1val1=ifelse( abs(A1_payoff)>=abs(A2_payoff), A1_payoff, A2_payoff), 
         opt1val2=ifelse( abs(A1_payoff)>=abs(A2_payoff), A2_payoff, A1_payoff),
         opt1prob1=ifelse(abs(A1_payoff)>=abs(A2_payoff), A1_prob, A2_prob), 
         opt1prob2=ifelse(abs(A1_payoff)>=abs(A2_payoff), A2_prob, A1_prob),
         opt2val1=ifelse( abs(B1_payoff)>=abs(B2_payoff), B1_payoff, B2_payoff), 
         opt2val2=ifelse( abs(B1_payoff)>=abs(B2_payoff), B2_payoff, B1_payoff),
         opt2prob1=ifelse(abs(B1_payoff)>=abs(B2_payoff), B1_prob, B2_prob), 
         opt2prob2=ifelse(abs(B1_payoff)>=abs(B2_payoff), B2_prob, B1_prob))

##* visualize variation in the data (values and probabilities) ----
p1=rieskamp_data %>%
  select(contains("val")) %>%
  gather() %>%
  ggplot(aes(x=value)) +
  geom_histogram(position="identity", bins=50) +theme_bw()+
  facet_wrap(.~key)

p2=rieskamp_data%>%
  select(contains("prob")) %>% select(starts_with("A") | starts_with("B")) %>% 
  gather() %>% 
  ggplot(aes(x=value)) +
  geom_histogram(position="identity", bins=50) + theme_bw()+
  facet_wrap(.~key)

ggarrange(p1, p2, ncol=2)

plot_hight <- 3
plot_width <- 7
save_plot <- purrr::partial(ggsave, width = plot_width, height =  plot_hight)
#save_plot(paste0(PATH_PLOTS, '/Rieskamp_val_probs_distr.pdf'))

##* build stan dataset ----
Nind=length(unique(rieskamp_data$subj_id))
stan.data=list(Nind=Nind,
               Nobs=nrow(rieskamp_data),
               id= pull(rieskamp_data, as.numeric(subj_id)),
               decision=pull(rieskamp_data, decision),
               opt1prob=pull(rieskamp_data, opt1prob1),
               opt2prob=pull(rieskamp_data, opt2prob1),
               opt1val1=pull(rieskamp_data, opt1val1),
               opt1val2=pull(rieskamp_data, opt1val2),
               opt2val1=pull(rieskamp_data, opt2val1),
               opt2val2=pull(rieskamp_data, opt2val2),
               nvar=5,
               u=matrix(rep(1, Nind), nrow=Nind), ndem=1)

## Run stan model ----
CPT_stanfit = stan(model_code = CPT_indcorr_two_outcome_gambles, data=stan.data, 
                     iter = 5000, warmup=3000, chains = 2, seed = 9000, verbose = TRUE, 
                     control = list(max_treedepth = 15, adapt_delta = 0.99))	

# save(CPT_stanfit, file=paste0(PATH_RESULTS, "/CPT_indcorr_rieskamp_two_outcome_gambles_data.RData"))

## Examine results----
##* elapsed time in hours----
rowSums(get_elapsed_time(CPT_stanfit))/3600
##* parameter estimates and traceplots----
if(IND){
  par_means=summary(CPT_stanfit, pars = c("mu_gamma","mu_delta", "mu_alpha", "mu_lambda", "mu_phi","lp__"), probs=c(0.025, 0.975))$summary
  ## reports the estimated correlation matrix 
  par_sd_cor_matrix=summary(CPT_stanfit, pars = c("sigma","Omega"), probs=c(0.025, 0.975))$summary
  traceplot(CPT_stanfit, pars=c("mu_gamma","mu_delta", "mu_alpha", "mu_lambda", "mu_phi","lp__"))
}else{
  par_means=summary(CPT_stanfit, pars = c("gamma","delta", "alpha", "lambda", "phi","lp__"), probs=c(0.025, 0.975))$summary
  traceplot(CPT_stanfit, pars=c("gamma","delta", "alpha", "lambda", "phi","lp__"))
}
par_means
par_sd_cor_matrix

#* model comparison and information criteria: log-pred density, WAIC----
ll= extract_log_lik(CPT_stanfit)
LPD_ind=waic(ll)[["estimates"]][1] + waic(ll)[["estimates"]][2]
WAIC_ind=waic(ll)

#* compare aggregate and individual models ----
model_comp=rbind(c("Aggregate model", in_sample_stats(stan.data, "agg_rieskamp_two_outcome_gambles_data" )),
                 c("Individual model", in_sample_stats(stan.data, "indcorr_rieskamp_two_outcome_gambles_data" )) )
colnames(model_comp)=c("Model" , "log-pred density", "WAIC", "Hit rate")
stargazer(model_comp, summary=F, align=F)   

## Posterior predictive checks ----
## needs both the aggregate and individual estimations (the stanfit objects) saved above
##* agg model ----
load(paste0(PATH_RESULTS, "/CPT_agg_rieskamp_two_outcome_gambles_data.RData") )
yrep_agg=t(as.matrix(CPT_stanfit,
                     pars = c("y_pred"))-1)
propRep_agg=tibble(prop_agg=rowMeans(yrep_agg))

gamma_agg  = summary(CPT_stanfit, pars = c("gamma"), probs=c(0.025, 0.975))$summary[,1]
delta_agg  = summary(CPT_stanfit, pars = c("delta"), probs=c(0.025, 0.975))$summary[,1]
alpha_agg  = summary(CPT_stanfit, pars = c("alpha"), probs=c(0.025, 0.975))$summary[,1]
lambda_agg = summary(CPT_stanfit, pars = c("lambda"), probs=c(0.025, 0.975))$summary[,1]

##* indcorr model ----
load(paste0(PATH_RESULTS, "/CPT_indcorr_rieskamp_two_outcome_gambles_data.RData") )
yrep_ind=t(as.matrix(CPT_stanfit,
                     pars = c("y_pred"))-1)
propRep_ind=tibble(prop_ind=rowMeans(yrep_ind))

gamma_ind  = summary(CPT_stanfit, pars = c("mu_gamma"), probs=c(0.025, 0.975))$summary[,1]
delta_ind  = summary(CPT_stanfit, pars = c("mu_delta"), probs=c(0.025, 0.975))$summary[,1]
alpha_ind  = summary(CPT_stanfit, pars = c("mu_alpha"), probs=c(0.025, 0.975))$summary[,1]
lambda_ind  = summary(CPT_stanfit, pars = c("mu_lambda"), probs=c(0.025, 0.975))$summary[,1]

##* correlation predicted and observed percentages of safer choices ----
rieskamp08_cors=rieskamp08_data %>%
  bind_cols(propRep_agg) %>%  bind_cols(propRep_ind) %>% 
  group_by(choicepair) %>% 
  mutate(prop_opt2=length(decision[decision==2])/n(),
         prop_opt2rep_ind=mean(prop_ind),
         prop_opt2rep_agg=mean(prop_agg)) %>% 
  ungroup() %>% 
  summarise(cor_agg=cor(prop_opt2, prop_opt2rep_agg, method="pearson"),
            cor_ind=cor(prop_opt2, prop_opt2rep_ind, method="pearson"))
rieskamp08_cors
 
##* plots of individual level parameters ----
params_est <- t(as.matrix(CPT_stanfit, pars = c("gamma", "delta", "alpha", "lambda", "phi")))
Niter=ncol(params_est)
Nind=30
nvar=5
colnames(params_est)=paste0(rep("niter", Niter), 1:Niter)

ind_params=as.data.frame(params_est)  %>% 
  mutate(param=c(rep("gamma", Nind), rep("delta", Nind),  
                 rep("alpha", Nind),  rep("lambda", Nind),
                 rep("phi", Nind)), 
         ID=rep(1:Nind, nvar) ) %>% 
  arrange(param, ID)
ind_params =ind_params %>% pivot_longer(starts_with("niter"), names_to="iter", values_to="value")

plot_ind_pars = ind_params %>%
  group_by(ID) %>% 
  mutate(median_value=median(value)) %>% 
  ungroup() %>% 
  ggplot(aes(x = reorder(ID, median_value), y = value))+ #color=param)) +
  stat_summary(fun = median,
               fun.min = function(z) {quantile(z,0.05)},
               fun.max = function(z) {quantile(z,0.95)},
               geom="pointrange", position = position_dodge(width=1))+
  labs(x = 'Participants', y = 'Parameter estimates', color="Parameter") +
  theme_bw()+theme(legend.position="bottom")+
  facet_wrap(~param)
plot_ind_pars


plot_hight <- 3
plot_width <- 8
save_plot <- purrr::partial(ggsave, width = plot_width, height =  plot_hight)
# save_plot(paste0(PATH_PLOTS, '/individual_median_estimates.pdf'))

##* Correlation between all individual-level parameters (medians) ----
df_corr= ind_params%>% 
  group_by(param, ID) %>% 
  summarise(medPar=median(value)) %>% 
  pivot_wider(names_from = param, values_from = medPar)
cor(as.matrix(df_corr[,-1]))

#* Plot estimated value and probability weighting functions ----
df_value=tibble(x= seq(-100, 100, 1) ) %>% 
  mutate(value_agg = ifelse(x>=0,  ( (abs(x))^alpha_agg ), -lambda_agg*( (abs(x))^alpha_agg )),
         value_ind = ifelse(x>=0,  ( (abs(x))^alpha_ind ), -lambda_ind*( (abs(x))^alpha_ind )), ) %>% 
  pivot_longer(-c(x), names_to="condition", values_to="Value") %>% 
  mutate(Model=if_else(grepl("agg", condition), "Aggregate",  "Individual") )

p1 = df_value %>%
  ggplot(aes(x=x, y=Value, color=Model)) +
  geom_line(size=0.5)+#smooth(se=F, method=lm, color="black", formula = y ~ splines::bs(x, 45), size=.5)+ 
  labs(x = "Value", y = "Utility", color="Model") +
  scale_color_manual(values=c("gray60", "darkgreen"))+
  theme_light()+theme(legend.position = "bottom")
p1

df_pwf=tibble(prob = seq(0, 1, 0.01) ) %>% 
  mutate(weight_agg = pwf_gonzalezwu(prob, delta_agg, gamma_agg),
         weight_ind = pwf_gonzalezwu(prob, delta_ind, gamma_ind) ) %>% 
  pivot_longer(-c(prob), names_to="condition", values_to="Weight") %>% 
  mutate(Model=if_else(grepl("agg", condition), "Aggregate",  "Individual") )

p2 = df_pwf %>%
  ggplot(aes(x=prob, y=Weight, color=Model)) +
  geom_line(size=0.5)+
  labs(x = "Probability", y = "Weight", color="Model") +
  ylim(c(-.001, 1.01))+
  scale_color_manual(values=c("gray60", "darkgreen"))+
  theme_light()+theme(legend.position = "bottom")
p2

plot_Value_PWF = grid.arrange(p1, p2, nrow = 1)
plot_Value_PWF

plot_height <- 4
plot_width <- 10
save_plot <- purrr::partial(ggsave, width = plot_width, height = plot_height)
# save_plot(paste0(PATH_PLOTS, "/Rieskamp_estimated_pwf_value_functions.pdf"), plot_Value_PWF)
