# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Read data and evaluate wind profiles to get roughness length (zo)
#  April 2024
#  This script was created by Anne Mendonça.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# free memory
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# set working directory
dir = "C:/PaperGRL"
setwd(dir)
# set data directory
data_dir = "Data/"   #input data directory
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# libraries
library(tidyverse)
theme_set(theme_bw())
library(scales)       # log ticks in axis
library(pals)
library(paletteer)
library(RColorBrewer)
library(lubridate)
library(ggpubr)       # gghistograms()
library(patchwork)    # plot_spacer()
library(ggforce)      # trans_reverser()
library(openair)
library(ggformula)
library(ggpmisc)
library(dplyr)
library(tidyr)
library(purrr)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# functions
# linear function to draw profile
linear_model = function(w,h,hc,p){
  d=p*hc # displacement height
  logh = log(h-d)
  f_model <- lm(logh ~ w)
  b=f_model$coefficients[1] #Intercept
  a=f_model$coefficients[2] #slope
  return(tibble(Intercept = b, slope = a))
}
#linear_model(dataf$windspeed,dataf$height)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# read data 
turb_stat <- read.csv("GRL/log_wind_profile.csv")
str(turb_stat)


#turb_stat <- na.omit(turb_stat) # remove NaN
turb_stat <- turb_stat[!is.na(turb_stat$wind_speed), ]

###################################################################

###  log(h-d)
hc=35 # canopy height
p=0.75 
# d=p*hc # displacement height
# turb_stat$height_log = log(turb_stat$height-d)


dataf = turb_stat %>%
  arrange(data_file, height) %>%
  group_by(data_file) %>% nest() %>% 
  mutate(model = map(data,~linear_model(.x$wind_speed,.x$height,hc,p)) ) %>%
  unnest(cols = c(data,model)) %>%
  mutate(z_o = exp(Intercept))
dataf <- as.data.frame(dataf)

dataf = dataf %>%
  mutate(mes = month(data_file))
sort(unique(dataf$mes))
#colnames(dataf)


# Select wind at 100 m
dataf = dataf %>%
  filter(height  %in% c(100) & mes %in% c(1,2,5,6,7,8,9,10,12) ) %>%
  dplyr::select(data_file,day_of_the_year,height,wind_dir_100,wind_speed,z_o,mes) %>% 
  rename("ws" = "wind_speed",
         "wd" = "wind_dir_100")

sort(unique(dataf$mes))

############## Figure Plot
### https://bookdown.org/david_carslaw/openair/sections/directional-analysis/polar-plots.html


filename = paste0("GRL/fig1a.png")
png(filename, width = 12, height = 9, bg = "white", units = "in",res=200)
polarPlot(
  dataf, 
  pollutant = "z_o",
  statistic = "mean",
  #statistic = "nwr",
  exclude.missing = TRUE,
  upper = 8, 
  limits = c(1, 3),
  k = 5, # 80
  key.header = expression("Mean " * z[o] * " (m)"),
  key.footer = "",
  key.position = "right",
  #par.settings=list(fontsize=list(text=30)),
  units = "m s-1",
  #breaks = seq(0, 10, by = 1)
)
dev.off()

