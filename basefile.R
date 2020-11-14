# 2020 November 14
# Greg
# 
# Phycocalidia tanegashimaensis
# タネガシマアマノリ
# 

##
library(tidyverse)
library(readxl)


# データはドロップボックスにある。
# Fig 1, 2, 5: 渡邉担当
# Fig. 3, 4, 6, 7, 8: GNN担当
dropbox = dir("~/Dropbox/MS_Watanabe_GNN_RT",
              pattern = "Phycocalidia",
              full.names = TRUE) 

data_folders = list.dirs(dropbox, full.names = TRUE, recursive = TRUE)
data_folders = data_folders[grepl("Fig_[1-8]", data_folders)]


################################################################################

path = dir(data_folders[grepl("Fig_3", data_folders)], full.names = TRUE,
           pattern = "xlsx")

f3data = tibble(fnames = path) %>% 
  mutate(data = map(fnames, read_xlsx, sheet = "Results"))

f3data = f3data %>% mutate(experiment = str_extract(fnames, "[0-9]{2}h")) %>% 
  select(experiment, data)

f3data = f3data %>% mutate(data = map(data, function(df) {
    df %>% rename_with(~c("Temp"), starts_with("...1"))
  }))


f3data = f3data %>% unnest() %>% select(experiment, Temp, starts_with("Y(II)"))

f3data = f3data %>% pivot_longer(cols = starts_with("Y(II)"),
                                 names_to = "number",
                                 values_to = "yield")

f3data = f3data %>% mutate(number = as.numeric(str_extract(number, "[0-9]+")))

f3data %>% group_nest(experiment) %>%
  mutate(out = walk(data, function(df) {print(tail(df))}))


# Check the data  and it looks ok.
 
ggplot(f3data) +
  geom_point(aes(x = Temp, y = yield),
             position = position_jitter(0.2),
             alpha = 0.5) +
  facet_grid(rows = vars(experiment))

write_csv(f3data, file = "f3data.csv")
