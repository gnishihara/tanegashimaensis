# 2020 November 14
# Greg
# 
# FIGURE 03
# Phycocalidia tanegashimaensis
# タネガシマアマノリ
# 

# Library
library(tidyverse)
library(readxl)
library(brms)
library(tidybayes)
library(bayesplot)

# Options
options(mc.cores = 8)

fname = "f3data.csv"
data = read_csv(fname)

data = data %>% rename(temperature = Temp)
data = data %>% mutate(ft = as_factor(temperature))

stan_function =  "
real ytmodel (real ps, real ha, real eta, real ktopt, real temperature) {
  real invkelvin = 1.0 / (temperature + 273.15);
  real gas_constant = 8.314/1000.0;
  real tmp0 = (1.0 / ktopt - invkelvin);
  real tmp1 = ha / gas_constant * tmp0 ;
  real tmp2 = (ha * eta) / gas_constant * tmp0;
  tmp1 = (ha * eta) * exp(tmp1);
  tmp2 = 1 - exp(tmp2);
  tmp2 = (ha * eta) - ha * tmp2;
  return ps * tmp1 / tmp2;
}
"

stanvars = stanvar(scode = stan_function, block = "functions")

# brms model 作成 ----
brmsmodel = brmsformula(yield ~ ytmodel(PS, HA, ET, KT, temperature),
                        PS ~ (1|experiment),
                        HA ~ (1|experiment),
                        ET ~ (1|experiment),
                        KT ~ (1|experiment),
                        zi ~ (1|ft),
                        nl = TRUE,
                        family = brms::zero_inflated_beta(link = "identity",
                                                          link_phi = "log",
                                                          link_zi = "logit"))

# Set the priors.

get_prior(brmsmodel, data = data)

PRIORS = 
  set_prior("normal(2, 10)",      class = "b", lb = 1, nlpar = "ET") +
  set_prior("normal(20+273, 50)", class = "b", lb = 273.15, nlpar = "KT") +
  set_prior("normal(0.5, 1)",     class = "b", lb = 0, nlpar = "PS") +
  set_prior("normal(80, 10)",     class = "b", lb = 0, nlpar = "HA")


CHAINS = 4
CORES  = CHAINS
SEED   = 2020
CONTROL = list(adapt_delta = 0.99, max_treedepth = 10)
testfit = brm(brmsmodel, 
                data = data,
                stanvars = stanvars, prior = PRIORS,
                iter = 10, chains = CHAINS, cores = CORES, seed = SEED)


fullfit = update(testfit, newdata = data,
                 iter = 2000,
                 control = CONTROL)

expose_functions(fullfit, vectorize = TRUE)

fullfit
# Prior predictive checks

y = data %>% pull(yield)
yrep = posterior_predict(fullfit, nsamples = 500)
ppc_dens_overlay(y,yrep)
 

newdata = fullfit$data %>% 
  expand(experiment, 
         temperature) %>% 
  mutate(ft = as_factor(temperature)) %>% 
  nest(data = everything()) %>% 
  mutate(expectation = map(data, 
                           add_linpred_draws, 
                           model = fullfit,
                           re_formula = NA)) %>% 
  mutate(expectation = map(expectation, mean_hdci)) %>% 
  select(-data) %>% 
  unnest(expectation)

ggplot() + 
  geom_line(aes(x = temperature, y = .value), data = newdata) +
  geom_ribbon(aes(x = temperature, ymin = .lower, ymax = .upper), data = newdata,
              alpha = 0.2) +
  geom_point(aes(x = temperature, y = yield, color = experiment),
             data = data,
             alpha = 0.5, position = position_jitter(0.2)) + 
  facet_grid(rows = vars(experiment))

