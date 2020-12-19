# 2020 November 14
# Greg
# 
# FIGURE 01
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

fname = "f1data.csv"
data = read_csv(fname)

stan_function =  "
real pemodel (real pm, real al, real rd, real light) {
  return pm * (1-exp(-al/pm * light)) - rd;
}
"

ggplot(data) +
  geom_point(aes(x = light, y = rate)) + 
  facet_wrap("temperature")


data %>% 
  group_by(temperature, light) %>% 
  summarise(vrate = sd(rate)) %>% 
  ggplot() +
  geom_point(aes(x = sqrt(light), y = vrate)) + 
  facet_wrap("temperature")

stanvars = stanvar(scode = stan_function, block = "functions")

# brms model 作成 ----
brmsmodel = brmsformula(rate ~ pemodel(PM, AL, RD, light),
                        PM ~ 1,
                        AL ~ 1,
                        RD ~ 1,
                        nl = TRUE)

# Set the priors.

get_prior(brmsmodel, data = data)

PRIORS = 
  set_prior("normal(50, 50)", class = "b", lb = 0, nlpar = "PM") +
  set_prior("normal(1, 10)",  class = "b", lb = 0, nlpar = "AL") +
  set_prior("normal(5, 10)",  class = "b", lb = 0, nlpar = "RD") 



CHAINS = 4
CORES  = CHAINS
SEED   = 2020
ITER = 2000
CONTROL = list(adapt_delta = 0.99, max_treedepth = 10)
SNAME = "figure01stanout.rds"

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
  tibble(temperature = name,
         model = c("f1", "f2", "f3"),
         bout = list(f1,f2,f3))
}

plan(list(
  tweak(multisession, workers = 4),
  tweak(multisession, workers = 4)
  ))


data = data %>% group_nest(temperature)

bout %<-%  run_models(data %>% pull(data), name = data$temperature)

bout %>% 
  mutate(ppc = map(bout, function(x) {
    y = x$data %>% pull(rate)
    yrep = posterior_predict(x, nsamples = 50)
    ppc_dens_overlay(y, yrep)
  })) %>% 
  pull(ppc) %>% 
  ggpubr::ggarrange(plotlist = .)

# Prior predictive checks

bout %>% slice(2) %>% pull(bout) 

bout = 
  bout %>% 
  mutate(fit = map(bout, function(x) {
    x$data %>% 
      complete(light = seq(min(light), max(light), length = 11)) %>% 
      add_fitted_draws(model = x) %>% 
      group_by(light) %>% 
      mean_hdci(.value)
  })) %>% 
  mutate(pred = map(bout, function(x) {
    x$data %>% 
      complete(light = seq(min(light), max(light), length = 11)) %>% 
      add_predicted_draws(model = x) %>% 
      group_by(light) %>% 
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



xlabel = "Irradiance (μmol photons m\\textsuperscript{-1} s\\textsuperscript{-1})"
ylabel = "Net photosynthesis  rates (μg O\\textsubscript{2} g\\rlap{\\textsuperscript{-1}}\\textsubscript{\\scriptsize{ww}} min\\textsuperscript{-1})"
xaxis = data %>% unnest(data) %>% pull(light) %>% range()
xaxis = range(pretty(xaxis))

p1 = 
ggplot() + 
  geom_line(aes(x = light, y = .value), 
            data = bout %>% unnest(fit)) +
  geom_ribbon(aes(x = light, ymin = .lower, ymax = .upper),
              data = bout %>% unnest(pred),
              alpha = 0.2) +
  geom_point(aes(x = light, y = rate),
             data = data %>% unnest(data),
             alpha = 0.5) + 
  geom_text(aes(x = Inf, y = Inf, label = sprintf("(%s)", l)),
            data = tibble(l = LETTERS[1:3],
                          temperature = data$temperature),
            hjust = 1, vjust = 1) +
  geom_text(aes(x = 1000, y = -25, label = sprintf("%s °C", temperature)),
            data = tibble(l = LETTERS[1:3],
                          temperature = data$temperature),
            hjust = 1, vjust = 0) +
  scale_x_continuous(xlabel, limits = xaxis, breaks = pretty(xaxis)) +
  scale_y_continuous(ylabel,
                     limits = c(-25, 100),
                     breaks = seq(-25, 100, by = 25)) + 
  facet_rep_grid(rows = vars(temperature)) +
  guides(color = F) +
  theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_blank())

wh = 80/25.4
ht = 1.5*80/25.4
pt = 10

# To run the tikz device, you will need to install tinytex and tikzxdevice,
# which is not trivial to accomplish.
filename = "./FIGURES/Figure01.tex"
tikz(filename, width = wh, height = ht, pointsize = pt, standAlone = TRUE)
print(p1)
dev.off()
tinytex::xelatex(filename)
img = magick::image_read(str_replace(filename, "tex", "pdf"), density = 600)
magick::image_write(img, str_replace(filename, "tex", "png"), format = "png")



