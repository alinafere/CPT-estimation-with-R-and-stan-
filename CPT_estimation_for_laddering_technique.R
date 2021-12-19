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
library(stargazer)

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
source(paste (PATH_FUNCTION,"/CPT_functions.R", sep="") ) 
##* compile stan model (individual or aggregate) ----
  if(IND){
    m = stan_model(model_code = CPT_indcorr_ladders)
  }else{  m = stan_model(model_code =CPT_agg_ladders)
  }  


#### Import Wu & Gonzalez 1996 data ----
wg96_data=read.csv(paste0(PATH_DATA,"/wu_gonzalez_1996_dataR.csv"), header = T)

### separate the value and probabilities in the lotteries
wg96_data_sep = wg96_data %>% 
  separate(r1, c("r1_prob","r1_value"), sep="([\\,\\$])",remove=F, extra="merge", convert = T) %>%
  separate(r2, c("r2_prob","r2_value"), sep="([\\,\\$])",remove=F, extra="merge", convert = T) %>%
  separate(s1, c("s1_prob","s1_value"), sep="([\\,\\$])",remove=F, extra="merge" , convert = T) %>%
  separate(r1_value, c("r1_vald","r1_val"), sep="([\\$])",remove=F, extra="merge", convert = T) %>% 
  separate(r2_value, c("r2_vald","r2_val"), sep="([\\$])",remove=F, extra="merge", convert = T) %>% 
  separate(s1_value, c("s1_vald","s1_val"), sep="([\\$])",remove=F, extra="merge", convert = T) %>% 
  mutate(s1_val=ifelse(is.na(s1_val), s1_vald, s1_val),
         s1_prob=ps )

### transform the choice variable into a decision variable
### choice: 1 - risky, 0 - safe
### decision: 1 - risky, 2 - safe;
wg96_data_sep = wg96_data_sep %>% 
  mutate(decision=ifelse(choice==0, 2, choice))

## create ID variable
wg96_data_sep=wg96_data_sep %>% group_by(subjectid) %>% 
  mutate( IND = cur_group_id() ) %>% ungroup() %>% 
  arrange(IND)

##* summary stats ----
nrow(wg96_data_sep) ## 4200 data points, with 5 ladders
wg96_data_sep %>% group_by(ladder) %>%  summarise(n()) ## 840 observations per ladder
length(unique(wg96_data_sep$IND)) ## 420 individuals
wg96_data_sep %>% group_by(IND) %>% ## 8 or 12 observations per individual across all ladders
  summarise(obs_per_individual = n()) %>% 
  ggplot(aes(x=obs_per_individual))+
  geom_histogram()

wg96_data_sep %>% group_by(IND, ladder) %>% ## 4 observations per individual within a ladder
  summarise(obs_per_individual = n()) %>% 
  ggplot(aes(x=obs_per_individual))+
  geom_histogram()

##* visualize data ----
wg96_data_sep %>% group_by(ladder, rung, s1_prob) %>%
  summarise(prop_safe=1-mean(choice)) %>%
  ggplot(aes(x=s1_prob, y=prop_safe))+
  geom_point()+
  scale_x_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))+
  labs(x = "Safer choice probability", y = "Safer choice proportion") +
  facet_wrap(~ladder, nrow=1, labeller = as_labeller( c("1"="Ladder 1", "2"="Ladder 2", "3"="Ladder 3", "4"="Ladder 4", "5"="Ladder 5")))+
  theme_light()+theme(legend.position="bottom")

plot_scale <- 3
plot_aspect <- 2.8
save_plot <- purrr::partial(ggsave, width = plot_aspect * plot_scale, height = 1 * plot_scale)
#save_plot(paste0(PATH_PLOTS, "/wg96_prop_safe_broken_by_Ladder.pdf"))

#* build stan dataset ----
## using full dataset here (but one could filter by ladder by outcommenting the filter() argument below)
stan_dataset=wg96_data_sep ## %>% filter(ladder %in% c(1, 4, 5))
Nind=length(unique(stan_dataset$IND))
stan_data=list(Nobs=nrow(stan_dataset),
              Nind=Nind,
              id=pull(stan_dataset, IND),
              decision=pull(stan_dataset, decision),
              opt1prob1=pull(stan_dataset, r1_prob),
              opt1prob2=pull(stan_dataset, r2_prob),
              opt2prob1=pull(stan_dataset, s1_prob),
              opt1val1= pull(stan_dataset, r1_val), 
              opt1val2= pull(stan_dataset, r2_val), 
              opt2val1= pull(stan_dataset, s1_val), 
              nvar=4, u=matrix(rep(1, Nind), nrow=Nind), ndem=1)
              
## Run stan model ----
CPT_stanfit = stan(model_code = CPT_indcorr_ladders, data=stan_data, 
                           iter = 5000, warmup=3000, chains = 2, seed = 9000, verbose = TRUE, 
                           control = list(max_treedepth = 15, adapt_delta = 0.99))	

## save stanfit object (outcomment below)
## save(CPT_stanfit, file=paste0(PATH_RESULTS, "/CPT_indcorr_wg96_ladders_data.RData"))

## Examine results----
##* elapsed time in hours----
rowSums(get_elapsed_time(CPT_stanfit))/3600
##* parameter estimates and traceplots----
if(IND){
  par_means=summary(CPT_stanfit, pars = c("mu_gamma","mu_delta", "mu_alpha", "mu_phi","lp__"), probs=c(0.025, 0.975))$summary
  ## reports the estimated correlation matrix 
  par_sd_cor_matrix=summary(CPT_stanfit, pars = c("sigma","Omega"), probs=c(0.025, 0.975))$summary
  traceplot(CPT_stanfit, pars=c("mu_gamma","mu_delta", "mu_alpha", "mu_phi","lp__"))
}else{
  par_means=summary(CPT_stanfit, pars = c("gamma","delta", "alpha", "phi","lp__"), probs=c(0.025, 0.975))$summary
  traceplot(CPT_stanfit, pars=c("gamma","delta", "alpha", "phi","lp__"))
}
par_means
par_sd_cor_matrix

#* model comparison and information criteria: log-pred density, WAIC----
ll= extract_log_lik(CPT_stanfit)
(LPD=waic(ll)[["estimates"]][1] + waic(ll)[["estimates"]][2])
(WAIC=waic(ll))

#* compare aggregate and individual models ----
model_comp=rbind(c("Aggregate model", in_sample_stats(stan_data, "agg_wg96_ladders_data" )),
                 c("Individual model", in_sample_stats(stan_data, "indcorr_wg96_ladders_data" )) )
colnames(model_comp)=c("Model" , "log-pred density", "WAIC", "Hit rate")
stargazer(model_comp, summary=F, align=F)   

## Posterior predictive checks ----
## needs both the aggregate and individual estimations (the stanfit objects) saved above
##* agg model ----
load(paste0(PATH_RESULTS, "/CPT_agg_wg96_ladders_data.RData"))
gamma_agg  = summary(CPT_stanfit, pars = c("gamma"), probs=c(0.025, 0.975))$summary[,1]
delta_agg  = summary(CPT_stanfit, pars = c("delta"), probs=c(0.025, 0.975))$summary[,1]
alpha_agg  = summary(CPT_stanfit, pars = c("alpha"), probs=c(0.025, 0.975))$summary[,1]

yrep_agg=t(as.matrix(CPT_stanfit,
                 pars = c("y_pred"))-1)
Niter=ncol(yrep_agg)
colnames(yrep_agg)=paste0(rep("yrep", Niter), 1:Niter)
yrep_df=as_tibble(yrep_agg) %>% 
  bind_cols(wg96_data_sep %>% 
              dplyr::select(s1_prob, ladder))

yrep_prop_agg = yrep_df %>% 
  pivot_longer(starts_with("yrep"), names_to="iter", values_to="value") %>% 
  group_by(s1_prob, ladder, iter) %>% 
  summarise(prop_safe_rep_agg = length(value[value==1])/n() )

##* indcorr model ----
load(paste0(PATH_RESULTS, "/CPT_indcorr_wg96_ladders_data.RData") )
gamma_ind  = summary(CPT_stanfit, pars = c("mu_gamma"), probs=c(0.025, 0.975))$summary[,1]
delta_ind  = summary(CPT_stanfit, pars = c("mu_delta"), probs=c(0.025, 0.975))$summary[,1]
alpha_ind  = summary(CPT_stanfit, pars = c("mu_alpha"), probs=c(0.025, 0.975))$summary[,1]

yrep_ind=t(as.matrix(CPT_stanfit,
                     pars = c("y_pred"))-1)
Niter=ncol(yrep_ind)
colnames(yrep_ind)=paste0(rep("yrep", Niter), 1:Niter)
yrep_df=as_tibble(yrep_ind) %>% 
  bind_cols(wg96_data_sep %>% 
              dplyr::select(s1_prob, ladder))

yrep_prop_ind = yrep_df %>% 
  pivot_longer(starts_with("yrep"), names_to="iter", values_to="value") %>% 
  group_by(s1_prob, ladder, iter) %>% 
  summarise(prop_safe_rep_ind = length(value[value==1])/n() )

##* plot pp checks ----
pp_checks_data =  wg96_data_sep %>% 
  group_by(s1_prob, ladder) %>%
  summarise(prop_safe=length(decision[decision==2])/n() ) %>% 
  left_join(yrep_prop_ind) %>% 
  left_join(yrep_prop_agg) %>% 
  pivot_longer(c(prop_safe_rep_ind, prop_safe_rep_agg), names_to="condition", values_to="prop") 

pp_checks_data %>% ggplot()+
  # ggplot(aes(x=s1_prob/100, y=prop_safe))+
  stat_summary(aes(x=s1_prob, y=prop), fun.data = mean_sdl, fun.args = list(mult = 1))+
  geom_point(aes(x=s1_prob, y=prop_safe), size=2, shape=18, position=position_dodge(width=.1))+
  geom_line(aes(x=s1_prob, y=prop_safe), linetype=3, color="gray60", position=position_dodge(width=.03))+
  scale_x_continuous(lim=c(0,1), breaks=c(0, 0.25, 0.5, 0.75, 1), labels=c("0", "0.25", "0.5", "0.75", "1"))+
  labs(x = "Safer choice probability", y = "Safer choice proportion") +
  facet_grid(condition~ladder, labeller = as_labeller( c("1"="Ladder 1", "2"="Ladder 2", "3"="Ladder 3", "4"="Ladder 4", "5"="Ladder 5",
                                                         "prop_safe_rep_ind"="Predicted (Indvidual)",  "prop_safe_rep_agg"="Predicted (Aggregate)") ) ) +
  theme_light()+theme(legend.position="bottom")

plot_scale <- 5
plot_aspect <- 3
save_plot <- purrr::partial(ggsave, width = plot_aspect * plot_scale, height = 1 * plot_scale)
# save_plot(paste0(PATH_PLOTS, '/pp_checks_props_safe_WG96_data.pdf'))

##* corr observed predicted props ----
pp_checks_corr=pp_checks_data %>% 
  group_by(condition, ladder, s1_prob, prop_safe) %>% 
  summarise(prop_mean=mean(prop)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=condition, values_from=prop_mean) %>% 
  summarise(cor_ind=round(cor(prop_safe, prop_safe_rep_ind, method="pearson"), 4),
            cor_agg=round(cor(prop_safe, prop_safe_rep_agg, method="pearson"), 4) ) 

stargazer(pp_checks_corr, summary=F, align=F, digits=4)

##* plots of individual level parameters ----
## Individual-level parameter estimates are unstable
## because there are very few data points per indiviudal (8 or 12 data points for 4 parameters)
params_est <- t(as.matrix(CPT_stanfit, pars = c("gamma", "delta", "alpha", "phi")))
Niter=ncol(params_est)
Nind=420
nvar=4
colnames(params_est)=paste0(rep("niter", Niter), 1:Niter)

ind_params=as.data.frame(params_est)  %>% 
  mutate(param=c(rep("gamma", Nind), rep("delta", Nind),  
                 rep("alpha", Nind), rep("phi", Nind)), 
         ID=rep(1:Nind, nvar) ) %>% 
  arrange(param, ID)
ind_params =ind_params %>% pivot_longer(starts_with("niter"), names_to="iter", values_to="value")

plot_ind_pars = ind_params %>% filter(ID %in% 1:100) %>% 
  filter( !(grepl("delta", param) | grepl("phi", param)) ) %>%
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

#* Plot estimated  value and probability weighting functions ----
df_value=tibble(x= seq(0, 100, 1) ) %>% 
  mutate(value_agg = ( (abs(x))^alpha_agg ),
         value_ind = ( (abs(x))^alpha_ind )) %>% 
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


