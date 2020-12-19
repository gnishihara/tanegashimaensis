# 2020 November 14
# Greg
# 
# FIGURE 02
# Phycocalidia tanegashimaensis
# タネガシマアマノリ
# 

# Library
library(tidyverse)
library(readxl)
library(rstan)
library(tidybayes)
library(bayesplot)
library(lemon)
library(ggpubr)
# Options
options(mc.cores = 8)

fname = "f2data.csv"
data = read_csv(fname)


ggplot(data) +
  geom_point(aes(x = temperature, y = rate)) + 
  facet_wrap("state")

data = data %>% 
  pivot_wider(names_from = state, values_from = rate) %>% 
  mutate(dark = -1 * dark) %>% arrange(temperature) %>% 
  filter(dark >0)

netp = data %>% pull(light)
resp = data %>% pull(dark)
temperature_idx = data %>% pull(temperature) %>% as.factor() %>% as.numeric()
temperature = data %>% pull(temperature) %>% unique()


#' The rest of the data to send to Stan.

nobs = length(netp) %>% as.integer()
ntemperature = length(temperature) %>% as.integer()
nmodeltemperature = 50 %>% as.integer()
modeltemperature = seq(min(temperature), 
                       max(temperature), 
                       length=nmodeltemperature)
Tcenter = 24

#' Setup the prior probability distribution
priormuEAmu = 60
priormuRDmu = 5
priormuPSmu = 10
priormuHAmu = 30
priormuETmu = 5
priormuKTmu = 12+273.15
priormuEAsigma = priormuEAmu
priormuRDsigma = priormuRDmu
priormuPSsigma = priormuPSmu
priormuHAsigma = priormuHAmu
priormuETsigma = priormuETmu
priormuKTsigma = priormuKTmu
priorCauchy    = c(1.0, 2.5)


CHAINS = 4
CORES = 4
ITER = 2000
SEED = 2020
CONTROL = list(adapt_delta = 0.95, max_treedepth = 11)

if(!file.exists("GPRP_Model.rds")) {
  smodel = stan_model(file = "GPRP_Model.stan", auto_write = TRUE)
} else {smodel = readRDS("GPRP_Model.rds")}

sout = sampling(smodel, chains = CHAINS, cores = CORES, iter = ITER, seed = SEED,
                control = CONTROL,
                show_messages = FALSE,
                verbose = FALSE,
                open_progress = FALSE,
                data = list(nobs = nobs, ntemperature = ntemperature,
                            temperature = temperature,
                            temperature_idx  =temperature_idx,
                            nmodeltemperature = nmodeltemperature,
                            modeltemperature = modeltemperature,
                            netp = netp, resp = resp,
                            Tcenter = Tcenter,
                            priorCauchy = priorCauchy,
                            priormuEAmu    = priormuEAmu,
                            priormuEAsigma = priormuEAsigma,
                            priormuPSmu    = priormuPSmu,
                            priormuPSsigma = priormuPSsigma,
                            priormuRDmu    = priormuRDmu,
                            priormuRDsigma = priormuRDsigma,
                            priormuHAmu    = priormuHAmu,
                            priormuHAsigma = priormuHAsigma,
                            priormuETmu    = priormuETmu,
                            priormuETsigma = priormuETsigma,
                            priormuKTmu    = priormuKTmu,
                            priormuKTsigma = priormuKTsigma))

predictedNP = tidy_draws(sout) %>% 
  select(contains("NPpred")) %>% 
  gather() %>% 
  group_nest(key) %>% 
  mutate(id = str_extract(key, "[0-9]+")) %>% 
  mutate(id = as.numeric(id)) %>% arrange(id) %>% 
  mutate(temperature = modeltemperature) %>% 
  unnest(data) %>% 
  group_by(temperature) %>% 
  mean_hdci(value)

expectedGP = tidy_draws(sout) %>% 
  select(contains("GPhat")) %>% 
  gather() %>% 
  group_nest(key) %>% 
  mutate(id = str_extract(key, "[0-9]+")) %>% 
  mutate(id = as.numeric(id)) %>% arrange(id) %>% 
  mutate(temperature = modeltemperature) %>% 
  unnest(data) %>% 
  group_by(temperature) %>% 
  mean_hdci(value)


predictedRP = tidy_draws(sout) %>% 
  select(contains("RPpred")) %>% 
  gather() %>% 
  group_nest(key) %>% 
  mutate(id = str_extract(key, "[0-9]+")) %>% 
  mutate(id = as.numeric(id)) %>% arrange(id) %>% 
  mutate(temperature = modeltemperature) %>% 
  unnest(data) %>% 
  group_by(temperature) %>% 
  mean_hdci(value)

################################################################################
# The plot.
# 
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
ylabel = "Photosynthesis and respiration rates (μg O\\textsubscript{2} g\\rlap{\\textsuperscript{-1}}\\textsubscript{\\scriptsize{ww}} min\\textsuperscript{-1})"

xaxis = temperature %>% range()
xaxis = range(pretty(xaxis))

p1 = ggplot() + 
  geom_line(aes(x = temperature, y = value), 
            data = predictedNP) +
  geom_ribbon(aes(x = temperature, 
                  ymin = .lower, ymax = .upper), 
              data = predictedNP,
              alpha = 0.2) +
  geom_point(aes(x = temperature, y = light),
             data = data,
             alpha = 0.5) + 
  annotate("text", x = 8, y = 0,
           label = "Net photosynthesis",
           vjust = 0, hjust = 0, size = 3) +
  scale_x_continuous("", 
                     limits = c(7.5, 36.5), 
                     breaks = temperature,
                     labels = rep("", length(temperature))) +
  scale_y_continuous("", limits = c(0, 80)) + 
  guides(color = F) +
  theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = NA))

p2 = ggplot() + 
  geom_line(aes(x = temperature, y = value), 
            data = expectedGP) +
  geom_ribbon(aes(x = temperature, 
                  ymin = .lower, ymax = .upper), 
              data = expectedGP,
              alpha = 0.2) +
  annotate("text", x = 8, y = 0,
           label = "Gross photosynthesis",
           vjust = 0, hjust = 0, size = 3) +
  scale_x_continuous("", 
                     limits = c(7.5, 36.5), 
                     breaks = temperature,
                     labels = rep("", length(temperature))) +
  scale_y_continuous(ylabel, limits = c(0, 80)) + 
  guides(color = F) +
  theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = NA))

p3 = ggplot() + 
  geom_line(aes(x = temperature, y = value), 
            data = predictedRP) +
  geom_ribbon(aes(x = temperature, 
                  ymin = .lower, ymax = .upper), 
              data = predictedRP,
              alpha = 0.2) +
  geom_point(aes(x = temperature, y = dark),
             data = data,
             alpha = 0.5) +
  annotate("text", x = 8, y = 20,
           label = "Dark respiration",
           vjust = 0, hjust = 0, size = 3) +
  scale_x_continuous(xlabel, 
                     limits = c(7.5, 36.5), 
                     breaks = temperature) +
  scale_y_reverse(   "", limits = c(20, 0)) + 
  guides(color = F) +
  theme_pubr()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.background = element_rect(fill = NA))

pfinal = ggarrange(p1,p2,p3, 
          ncol = 1,
          align = "hv",
          labels = c("(A)", 
                     "(B)", 
                     "(C)"),
          label.x = 1, 
          label.y = 1,
          hjust = 1,
          vjust = 1,
          font.label = list(face = "plain"))


wh = 80/25.4
ht = 1.5*80/25.4
pt = 10

# To run the tikz device, you will need to install tinytex and tikzxdevice,
# which is not trivial to accomplish.
filename = "./FIGURES/Figure02.tex"
tikz(filename, width = wh, height = ht, pointsize = pt, standAlone = TRUE)
print(pfinal)
dev.off()
tinytex::xelatex(filename)

img = magick::image_read(str_replace(filename, "tex", "pdf"), density = 600)
img = img %>% magick::image_background(color = "white")
magick::image_write(img, str_replace(filename, "tex", "png"), format = "png")



