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
library(lemon)
library(ggpubr)
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

ggplot(data) +
  geom_point(aes(x = temperature, y = yield)) + 
  facet_wrap("experiment")

stanvars = stanvar(scode = stan_function, block = "functions")

# brms model 作成 ----
brmsmodel = brmsformula(yield ~ ytmodel(PS, HA, ET, KT, temperature),
                        PS ~ 1,
                        HA ~ 1,
                        ET ~ 1,
                        KT ~ 1,
                        zi ~ temperature,
                        nl = TRUE,
                        family = brms::zero_inflated_beta(link = "identity",
                                                          link_phi = "log",
                                                          link_zi = "logit"))

# Set the priors.

# get_prior(brmsmodel, data = data)

PRIORS = 
  set_prior("normal(2, 10)",      class = "b", lb = 1, nlpar = "ET") +
  set_prior("normal(20+273, 50)", class = "b", lb = 273.15, nlpar = "KT") +
  set_prior("normal(0.8, 2)",     class = "b", lb = 0, nlpar = "PS") +
  set_prior("normal(80, 10)",     class = "b", lb = 0, nlpar = "HA")


CHAINS = 4
CORES  = CHAINS
SEED   = 2020
ITER = 2000
CONTROL = list(adapt_delta = 0.99, max_treedepth = 10)
SNAME = "figure03stanout.rds"

# if(file.exists(SNAME)) {file.remove(SNAME)}

if(!file.exists(SNAME)){
  testfit = brm(brmsmodel, 
                data = data,
                stanvars = stanvars, prior = PRIORS,
                iter = 10, chains = CHAINS, cores = CORES, seed = SEED,
                silent = TRUE,
                refresh = 0,
                open_progress = FALSE)
  
  save(testfit, file = SNAME)
} else {
  load(SNAME)
}

expose_functions(testfit, vectorize = TRUE)
library(future)

run_models = function(data, name) {
  ITER = 2000
  CONTROL = list(adapt_delta = 0.99, max_treedepth = 10)
  f1 = update(testfit, newdata = data[[1]], iter = ITER, seed = 2020, control = CONTROL, future = TRUE, silent=TRUE, refresh=0,open_progress=FALSE)
  f2 = update(testfit, newdata = data[[2]], iter = ITER, seed = 2020, control = CONTROL, future = TRUE, silent=TRUE, refresh=0,open_progress=FALSE)
  f3 = update(testfit, newdata = data[[3]], iter = ITER, seed = 2020, control = CONTROL, future = TRUE, silent=TRUE, refresh=0,open_progress=FALSE)
  f4 = update(testfit, newdata = data[[4]], iter = ITER, seed = 2020, control = CONTROL, future = TRUE, silent=TRUE, refresh=0,open_progress=FALSE)
  tibble(experiment = name,
         model = c("f1", "f2", "f3", "f4"),
         bout = list(f1,f2,f3,f4))
}

plan(list(
  tweak(multisession, workers = 4),
  tweak(multisession, workers = 4)
  ))


data =data %>% group_nest(experiment)

bout %<-%  run_models(data %>% pull(data), data$experiment)

bout %>% 
  mutate(ppc = map(bout, function(x) {
    y = x$data %>% pull(yield)
    yrep = posterior_predict(x, nsamples = 50)
    ppc_dens_overlay(y, yrep)
  })) %>% 
  pull(ppc) %>% 
  ggpubr::ggarrange(plotlist = .)

# Prior predictive checks


bout = bout %>% 
  mutate(fit = map(bout, function(x) {
    x$data %>% 
      complete(temperature = seq(min(temperature), max(temperature))) %>% 
      add_fitted_draws(model = x) %>% 
      group_by(temperature) %>% 
      mean_hdci(.value)
  })) %>% 
  mutate(pred = map(bout, function(x) {
    x$data %>% 
      complete(temperature = seq(min(temperature), max(temperature))) %>% 
      add_predicted_draws(model = x) %>% 
      group_by(temperature) %>% 
      mean_hdci(.prediction)
  }))


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



xlabel = "Temperature (°C)"
ylabel = "Effective quantum yield"
xaxis = data %>% unnest(data) %>% pull(temperature) %>% range()
xaxis = range(pretty(xaxis))
p1 = ggplot() + 
  geom_line(aes(x = temperature, y = .value), 
            data = bout %>% unnest(fit)) +
  geom_ribbon(aes(x = temperature, ymin = .lower, ymax = .upper), 
              data = bout %>% unnest(fit),
              alpha = 0.6) +
  # geom_ribbon(aes(x = temperature, ymin = .lower, ymax = .upper), 
  #             data = bout %>% unnest(pred),
  #             alpha = 0.2) +
  geom_point(aes(x = temperature, y = yield),
             data = data %>% unnest(data),
             alpha = 0.5, position = position_jitter(0.2)) + 
  geom_text(aes(x = Inf, y = Inf, label = sprintf("(%s)", l)),
            data = tibble(l = LETTERS[1:4],
                          experiment = data$experiment),
            hjust = 1, vjust = 1) +
  geom_text(aes(x = 0, y = 0, label = sprintf("%s h", experiment2)),
            data = tibble(l = LETTERS[1:4],
                          experiment = data$experiment) %>% 
              mutate(experiment2 = str_extract(experiment, "[0-9]+")),
            hjust = 0, vjust = 0) +
  scale_x_continuous(xlabel, limits = xaxis, breaks = pretty(xaxis)) +
  scale_y_continuous(ylabel, limits = c(0, 0.8)) + 
  facet_rep_grid(rows = vars(experiment)) +
  guides(color = F) +
  theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_blank())

wh = 80/25.4
ht = 80/25.4
pt = 10

# To run the tikz device, you will need to install tinytex and tikzxdevice,
# which is not trivial to accomplish.
filename = "./FIGURES/Figure03.tex"
tikz(filename, width = wh, height = 2*ht, pointsize = pt, standAlone = TRUE)
print(p1)
dev.off()
tinytex::xelatex(filename)

img = magick::image_read(str_replace(filename, "tex", "pdf"), density = 600)

magick::image_write(img, str_replace(filename, "tex", "png"), format = "png")



