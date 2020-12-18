# 2020 November 14
# Greg
# 
# FIGURE 04
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

fname = "f4data.csv"
data = read_csv(fname)

data = data %>% rename(temperature = Temperature,
                       value = sr)
data = data %>% mutate(ft = as_factor(temperature)) %>% 
  filter(is.finite(value))

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
brmsmodel = brmsformula(value ~ ytmodel(PS, HA, ET, KT, temperature) + Z,
                        PS ~ 1,
                        HA ~ 1,
                        ET ~ 1,
                        KT ~ 1,
                        Z ~ 1,
                        nl = TRUE,
                        family = gaussian)

# Set the priors.

get_prior(brmsmodel, data = data)

PRIORS = 
  set_prior("normal(2, 10)",      class = "b", lb = 1, nlpar = "ET") +
  set_prior("normal(20+273, 50)", class = "b", lb = 273.15, nlpar = "KT") +
  set_prior("normal(0.5, 1)",     class = "b", lb = 0, nlpar = "PS") +
  set_prior("normal(80, 10)",     class = "b", lb = 0, nlpar = "HA") +
  set_prior("normal(0, 2)", class = "b", nlpar = "Z")


CHAINS = 4
CORES  = CHAINS
SEED   = 2020
CONTROL = list(adapt_delta = 0.99, max_treedepth = 10)
SNAME = "figure04stanout.rds"

if(!file.exists(SNAME)){
  testfit = brm(brmsmodel, 
                data = data,
                stanvars = stanvars, prior = PRIORS,
                iter = 10, chains = CHAINS, cores = CORES, seed = SEED)
  
  
  fullfit = update(testfit, newdata = data,
                   iter = 2000,
                   seed = SEED,
                   control = CONTROL)
  
  save(fullfit, file = SNAME)
} else {
  load(SNAME)
}

expose_functions(fullfit, vectorize = TRUE)

fullfit
# Prior predictive checks

y = data %>% pull(value)
yrep = posterior_predict(fullfit, nsamples = 500)
ppc_dens_overlay(y,yrep)


# Publication ready figure
library(ggpubr)
library(tikzDevice)
library(tinytex)
options(tikzXelatexPackages = c(
  "\\usepackage{xeCJK}\n",
  "\\usepackage{tikz}\n",
  "\\usepackage{fontawesome}\n",
  "\\usepackage[active,tightpage,xetex]{preview}\n",
  "\\usepackage{unicode-math,xunicode}\n",
  "\\PreviewEnvironment{pgfpicture}\n",
  "\\setlength\\PreviewBorder{0pt}\n",
  "\\defaultfontfeatures{Ligatures=TeX,Scale=MatchLowercase}\n",
  "\\setmainfont[]{Noto Serif}\n",
  "\\setsansfont[]{Noto Serif}\n",
  "\\setmonofont[Mapping=tex-ansi]{Source Code Pro}\n",
  "\\setCJKmainfont[]{Noto Serif CJK JP}\n",
  "\\setmathfont{XITS Math}\n",
  "\\usepackage{microtype}\n",
  "\\UseMicrotypeSet[protrusion]{basicmath}\n",
  "\\usepackage{textgreek}\n",
  "\\usetikzlibrary{matrix, tikzmark, fit, shapes}\n",
  "\\usetikzlibrary{arrows, arrows.meta, shapes.geometric, positioning}\n",
  "\\usetikzlibrary{calc, intersections,through}\n",
  "\\usetikzlibrary{decorations.text, decorations.markings}\n"
),
tikzDefaultEngine="xetex",
tikzMetricPackages = c("\\usetikzpackage{calc}"))


data = fullfit$data
newdata = fullfit$data %>% 
  expand(temperature =seq(min(temperature), max(temperature), by = 0.25)) %>% 
  mutate(ft = as_factor(temperature)) %>% 
  nest(data = everything()) %>% 
  mutate(expectation = map(data, 
                           add_linpred_draws, 
                           model = fullfit,
                           re_formula = NA)) %>% 
  mutate(expectation = map(expectation, mean_hdci)) %>% 
  select(-data) %>% 
  unnest(expectation)

xlabel = "Temperature (°C)"
ylabel = "Relative growth rate"
xlimits = data %>% pull(temperature) %>% unique()

p1 = 
  ggplot() + 
  geom_line(aes(x = temperature, y = .value), data = newdata) +
  geom_ribbon(aes(x = temperature, ymin = .lower, ymax = .upper), data = newdata,
              alpha = 0.2) +
  geom_point(aes(x = temperature, y = value), data = data, alpha = 0.5) +
  scale_x_continuous(xlabel, limits = range(xlimits),
                     breaks = xlimits) +
  scale_y_continuous(ylabel, limits = c(-0.1, 0.3)) +
  theme_pubr() +
  theme(legend.background = element_blank(),
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.direction = "horizontal",
        legend.title = element_blank())
p1
wh = 80/25.4
ht = 80/25.4
pt = 10

# To run the tikz device, you will need to install tinytex and tikzxdevice,
# which is not trivial to accomplish.
filename = "./FIGURES/Figure04.tex"
tikz(filename, width = wh, height = ht, pointsize = pt, standAlone = TRUE)
print(p1)
dev.off()
tinytex::xelatex(filename)

img = magick::image_read(str_replace(filename, "tex", "pdf"), density = 600)

magick::image_write(img, str_replace(filename, "tex", "png"), format = "png")


datatable = fullfit %>% fixef(summary = F) %>% 
  as_tibble() %>% 
  mutate(HD_Intercept = HA_Intercept * ET_Intercept,
         CT_Intercept = KT_Intercept - 273.15) %>% 
  select(CT = CT_Intercept,
         HA = HA_Intercept,
         HD = HD_Intercept,
         PS = PS_Intercept,
         Z = Z_Intercept) %>% 
  gather() %>% 
  group_by(key) %>% 
  mean_hdci()

