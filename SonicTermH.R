# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  Read data and evaluate means and standard deviations
#  April 2024
#  This script was created by Anne Mendonça with the support of Luca Mortarini.
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
library(Hmisc)   # "A" %nin% "B"
library(data.table)
library(dplyr)
library(tidyr)
library(tidync)
library(ggpmisc)
library(egg)
library(minpack.lm)
library(xtable)

##### functions  ###############################################################
# function to find zi by flux ###########
getCLheight <- function(flux, height) {
  zi_CL <- 0
  zi_flux <- 0
  zi_CL1 <- 0
  zi_CL2 <- 0
  
  for (j in 2:(length(flux))) { # when the flux reaches 5%
    if (!is.na(flux[j]) && zi_CL1 == 0 && flux[j] <= 0.06){
      zi_CL1 <- height[j] 
      
    }
  }
  
  for (j in 2:(length(flux)-1)) { # Doesn't reach 5%
    dif = abs(flux[j] - flux[j+1])
    if (!is.na(dif) && zi_CL2 == 0 && dif <= 0.005){
      zi_CL2 <- height[j] 
      
    }
  }
  
  
  zi_CL <- if (zi_CL1 != 0 && zi_CL2 != 0) {
    min(zi_CL1, zi_CL2)
  } else if (zi_CL1 != 0) { zi_CL1
  } else if (zi_CL2 != 0) { zi_CL2
  } else {   0
  }
  
  if (zi_CL == 0) { # If zi > tower height (316m)
    zi_flux <- 0 
    zi_CL <- 316
  } else {
    zi_flux <- flux[which(height==zi_CL)]
  }
  
  return(data.frame(zi_flux = zi_flux, zi_CL = zi_CL))
}

getZiLenschow <- function(flux, height) {
  zi_CL <- 0
  
  flux_data = arrange(tibble(height = height,flux = flux),height)
  
  for (j in 2:(length(flux))) { # Start loop to find the height of CL
    if (zi_CL == 0 & flux_data$flux[j] <= 0.06) {
      zi_CL <- flux_data$height[j] }
  }
  
  if (zi_CL == 0) { # Se nao encontrar fluxo <= 0.06
    zi_CL <- NA
  }
  
  return(zi_CL)
}

fitflux_Lenschow <- function(flux, height, h_c, zi, name) { # Estima h
  
  if (name == "uw") r = 1.75
  if (name == "wT") r = 1.5
  
  index_zi = which(height <= zi)
  flux = flux[index_zi]
  height = height[index_zi]
  
  index_flux = which(flux > 0)
  flux = flux[index_flux]
  height = height[index_flux]
  flux_data = arrange(tibble(height = height,flux = flux),height)
  index_50 = which(flux_data$height == 50)
  
  if(flux_data$flux[index_50] > 1) {
    h_c = 50
    index_h = which(flux_data$height == max(flux_data$height))
    height = flux_data$height[index_50:index_h]
    flux   = flux_data$flux[index_50:index_h]/flux_data$flux[index_50]
    #    print(c(index_50,index_h))
  }
  
  
  model <- function(params, h) { # Modelo matemático
    h_value <- params[1] # parâmetro h, a ser ajustado
    (1 - (h - h_c) / h_value)^r
  }
  
  params <- c(zi)
  
  fn <- function(params, h, flux) { # Calcula os resíduos
    model_output <- model(params, h)
    flux - model_output
  }
  
  # nls.lm faz o ajuste não linear usando o algoritmo de Levenberg-Marquardt
  result <- nls.lm(par = params, fn = fn, h = height, flux = flux)
  h_found = as.numeric(result$par[1])
  #  
  #    hfit = seq(min(height), zi, 1)
  #    FIT = (1-(hfit-h_c)/h_found)**r_found
  return(tibble(h_teo = h_found, r_teo = r))
}

fitflux_LM <- function(flux, height, h_c, zi, name) { # Estima h e r
  
  if (name == "uw") rstart = 1.75
  if (name == "wT") rstart = 1.5
  
  index_zi = which(height <= zi)
  flux = flux[index_zi]
  height = height[index_zi]
  
  index_flux = which(flux > 0)
  flux = flux[index_flux]
  height = height[index_flux]
  flux_data = arrange(tibble(height = height,flux = flux),height)
  index_50 = which(flux_data$height == 50)
  
  if(flux_data$flux[index_50] > 1) {
    h_c = 50
    index_h = which(flux_data$height == max(flux_data$height))
    height = flux_data$height[index_50:index_h]
    flux   = flux_data$flux[index_50:index_h]/flux_data$flux[index_50]
    #    print(c(index_50,index_h))
  }
  
  
  model <- function(params, h) { # Modelo matemático
    h_value <- params[1] # parâmetro h, a ser ajustado
    r_value <- params[2] # parâmetro r, a ser ajustado
    (1 - (h - h_c) / h_value)^r_value
  }
  
  params <- c(zi, rstart)
  
  fn <- function(params, h, flux) { # Calcula os resíduos
    model_output <- model(params, h)
    flux - model_output
  }
  
  # nls.lm faz o ajuste não linear usando o algoritmo de Levenberg-Marquardt
  result <- nls.lm(par = params, fn = fn, h = height, flux = flux)
  h_found = as.numeric(result$par[1])
  r_found = as.numeric(result$par[2])
  #  
  #    hfit = seq(min(height), zi, 1)
  #    FIT = (1-(hfit-h_c)/h_found)**r_found
  return(tibble(h = h_found, r = r_found, h_c = h_c))
}

correct_norm = function(flux,height){
  index_50 = which(height == 50)
  
  if(flux[index_50] > 1) {
    flux   = flux/flux[index_50]
  }
  
  return(correct_flux = flux)
}
##### parameters ###############################################################
h_canopy = 35
threshold <- 1.5
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Set color palette
ATTO_palette = c(stepped(20)[9], stepped3(20)[1], stepped(20)[2],stepped3(20)[5],stepped2(20)[3] )
pal.bands(ATTO_palette)

####################### Read sonic & termohig data ###############################

sonic_termh <- read.csv("GRL/data_sonic_termh_2022.csv")

# Potential temperature: "mean_Tp"
sonic_termh$mean_Tp = sonic_termh$mean_Th + ((9.8/1005)*sonic_termh$height)

# Potential temperature difference from the value observed at 298 m
turb_rfp_298 = sonic_termh %>%
  ungroup() %>%
  filter(height == 298) %>%
  dplyr::select(data_file, mean_Tp) %>%
  rename(mean_Tp298 = mean_Tp)

# join (merge) 
sonic_termh <- sonic_termh %>%
  left_join(turb_rfp_298, by = "data_file") %>%
  mutate(dif_Tp298 = mean_Tp - mean_Tp298) %>% ## Diference T(x)-T(298)
  select(-mean_Tp298)  
rm(turb_rfp_298)

############################################################

# Wind direction classes: "wind_pred"
data_WD <- sonic_termh %>%
  mutate(wind_pred = case_when(
    wind_dir_100 >= 0 & wind_dir_100 <= 70 ~ "a.00-70",
    wind_dir_100 >= 90 & wind_dir_100 <= 180 ~ "b.90-180",
    TRUE ~ "Others"  
  ))
table(factor(data_WD$wind_pred))
range(data_WD$wind_dir_100[data_WD$wind_pred %in% c("a.00-70")])
range(data_WD$wind_dir_100[data_WD$wind_pred %in% c("b.90-180")])

# Stability classes: "stability"
data_WD <- data_WD %>%
  mutate(stability = case_when(
    ZonL_35 >= 0 & ZonL_35 < 0.2 ~ "a.NN",
    ZonL_35 >= 0.2 & ZonL_35 < 0.4 ~ "b.NN-WS",
    ZonL_35 >= 0.4 & ZonL_35 < 0.6 ~ "c.WS-1",
    ZonL_35 >= 0.6 & ZonL_35 < 1.0 ~ "d.WS-2",
    ZonL_35 >= 1.0 & ZonL_35 < 2.0 ~ "e.WS-VS",
    ZonL_35 >= 2.0 & ZonL_35 < 10.0 ~ "f.VS",
    TRUE ~ "Others" 
  ))
table(factor(data_WD$stability))
range(data_WD$ZonL_35[data_WD$stability %in% c("a.NN")])
range(data_WD$ZonL_35[data_WD$stability %in% c("c.WS-1")])

############### remove outlier by wT ################

datanew_wT <- data_WD %>% 
  group_by(stability, height, wind_pred) %>%
  mutate(lower_limit_wT = quantile(wT, 0.25) - threshold * IQR(wT),
         upper_limit_wT = quantile(wT, 0.75) + threshold * IQR(wT)) %>%
  filter (wT >= lower_limit_wT & wT <= upper_limit_wT) %>%
  ungroup()

################# Average by classes (stability & wind_pred) #########

datanew_wT_stat <- datanew_wT %>% 
  group_by(height, stability, wind_pred) %>%
  summarise(mean_wT = mean(wT_norm, na.rm = TRUE), median_wT = median(wT_norm, na.rm = TRUE),
            mean_uw = mean(uw_norm, na.rm = TRUE), median_uw = median(uw_norm, na.rm = TRUE),
            mean_ws = mean(wind_speed, na.rm = TRUE), median_ws = median(wind_speed, na.rm = TRUE),
            mean_Th = mean(mean_Th, na.rm = TRUE), median_Th = median(mean_Th, na.rm = TRUE),
            mean_Tp = mean(mean_Tp, na.rm = TRUE), median_Tp = median(mean_Tp, na.rm = TRUE),
            mean_dif_Tp = mean(dif_Tp298, na.rm = TRUE), median_dif_Tp = median(dif_Tp298, na.rm = TRUE),
            .groups = 'keep'
  )


############### facet wT and uw by wind direction: FIG 1b-g ##################

flux_values = datanew_wT_stat %>%
  filter (wind_pred != "Others") %>%
  group_by(height, stability, wind_pred, median_wT, median_uw ) %>%
  pivot_longer(cols = c(median_wT, median_uw),
               names_to = "variable", values_to = "value") %>%
  select(height, stability, wind_pred, variable, value ) %>%
  mutate(z_on_h = height/35,
         value = ifelse(height == 316, NA, value))

names = c("median_wT"=expression(paste(bar("w'T'"))), 
          "median_uw"=expression(paste(bar("u'w'")))   )
colr = c("median_wT"="#E41A1C", "median_uw"="#2166AC")
liness = c("a.00-70" = "solid", "b.90-180" = "dashed")
winddir = c("a.00-70" = "lower", "b.90-180" = "higher")

fig1 = 
  ggplot(flux_values, aes(x = height, y = value)) + 
  geom_point(size = 2, aes(color = variable)) +
  geom_line(aes(group = interaction(wind_pred, variable),
                linetype = wind_pred, color = variable), linewidth = .5) +  
  geom_vline(xintercept = 35, lty = 2) + 
  geom_hline(yintercept = 0, lty = 2) +
  xlab(expression("Height (m)")) +
  ylab(expression(paste(bar("w'T'") / bar("w'T'")[35]," ,   ",bar("u'w'") / bar("u'w'")[35]  ))) +
  scale_x_continuous(breaks = c(35,100,200,300)) +
  scale_color_manual(values = colr, name = "Profile", labels = names) +
  scale_linetype_manual(values = liness, name = "Roughness", labels = winddir) +
  scale_y_continuous(limits = c(-.5, 1.5)) +
  coord_flip()+ theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 16),
        strip.text = element_text(face="bold",size = 14),
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14),
        #legend.position = c(0.92, 0.85),  # Posição da legenda (x, y)
        legend.background = element_rect( colour = "grey"),
        strip.background = element_rect(colour="white", fill="white")) +
  facet_wrap(. ~ stability, ncol = 3, scales = "fixed") 
fig1

ggsave(plot = (fig1) ,
       filename = "GRL/fig1_b-g.png", 
       width = 11, height = 7, units = "in", scale = 1, dpi = 300) 

rm(flux_values, names, colr, liness, winddir, fig1)


############### Tp and ws by wind direction: FIG 2 ##################

ws_Tp = datanew_wT_stat %>%
  filter (wind_pred != "Others") %>%
  group_by(wind_pred, stability, height, mean_dif_Tp, mean_ws) %>%
  pivot_longer(cols = c(mean_dif_Tp, mean_ws),
               names_to = "variable", values_to = "value") %>%
  select(wind_pred, stability, height, variable, value )%>%
  mutate(z_on_h = height/35,
         value = ifelse(height == 316, NA, value))

# Removing NAN
ws_Tp <- ws_Tp[complete.cases(ws_Tp[, c("value")]), ]

# Plot
winddir = c("a.00-70" = "Lower roughness", "b.90-180" = "Higher roughness")
colr = c("#B2182B","#FD8D3C","#0CB702","#00A9FF","#C77CFF","#FF68A1")
var = c("mean_dif_Tp" = "Potential temperature", "mean_ws" = "Wind speed")
lb = c("NN", "NN-WS", "WS-1", "WS-2", "WS-VS", "VS")

fig2 =
  ggplot(ws_Tp, 
         aes(x = height, y = value, group = factor(stability))) + 
  geom_point(size = 2, aes(color = factor(stability))) +
  geom_line(aes(group = factor(stability), color = factor(stability)), linewidth = .5) +  
  geom_hline(yintercept = 0, lty = 2) +
  xlab(expression("Height (m)")) +
  ylab(expression("         "~ T[(z)] - T[(298)] ~ "        " ~ "                Absolute value (m s"^{-1}~")")) +
  scale_color_manual(name = "Stability", values = colr, labels = lb) +
  scale_x_continuous(breaks = c(35,100,200,300)) +
  coord_flip()+ theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, color = "black"), # Rótulos dos eixos x e y em preto
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        strip.text = element_text(face="bold",size = 14),
        legend.text = element_text(size = 14), legend.title = element_text(size = 14),
        legend.background = element_rect( colour = "grey"),
        strip.background = element_rect(colour="white", fill="white")) +
  facet_grid( wind_pred ~ variable, scales = "free" ,
              labeller = labeller(variable = as_labeller(var),
                                  wind_pred = as_labeller(winddir))) 
fig2

# Add text a), b), ..., to each panel of the figure
labels <- data.frame(
  wind_pred = rep(names(winddir), length(var)),
  variable = rep(names(var), each = length(winddir)),
  label = paste0(letters[1:(length(winddir) * length(var))], ")" ),
  stability = rep(unique(sort(ws_Tp$stability)), length(var))
)

fig2_ = fig2 + 
  geom_text(data = labels, aes(label = label), 
            x = Inf, y = Inf, hjust = 14, vjust = 1.5, size=5)
fig2_

ggsave(plot = (fig2_)   +
         plot_layout(guides = "collect"),
       filename = "GRL/fig2.png",
       width = 8, height = 7, units = "in", scale = 1, dpi=300) 
rm(fig2, fig2_, ws_Tp, winddir, var, colr, lb, labels)


################### Get zi by average profiles ################

# Loop by classes: stability x wind_pred
df_raw <- data.frame() 

for (mn in sort(unique(datanew_wT_stat$wind_pred[datanew_wT_stat$wind_pred != "Others"]))) { 
  # mn = sort(unique(datanew_wT_stat$wind_pred))[1]
  dataz = datanew_wT_stat %>% filter (wind_pred == mn)
  print(paste0("Loop by wind dir_",dataz$wind_pred[1]))
  
  for (st in sort(unique(dataz$stability))) { 
    # st = sort(unique(dataz$stability))[1]
    datay = dataz %>% filter (stability == st)
    print(datay$stability[1])
    
    datay <- datay[order(datay$height),] 
    
    ### Estimar zi pelo perfil mediano - Rawdata
    result_wT <- getCLheight (datay$median_wT, datay$height) 
    result_uw <- getCLheight (datay$median_uw, datay$height)
    
    # Salvar estimativas de zi
    datay$zi_wT = result_wT$zi_CL
    datay$zi_uw = result_uw$zi_CL
    
    ### Average profiles with zi estimate values
    df_raw <- rbind(df_raw, datay)
    
  }
}

rm(datay,result_wT,result_uw,dataz,mn,st)

# Check zi estimates
data_zi <- df_raw %>%
  filter(height != 316)%>%
  dplyr::select(c(height, wind_pred, stability, median_uw, median_wT, zi_wT, zi_uw)) %>%
  rename(uw = median_uw,
         wT = median_wT) %>%
  mutate(height = as.numeric(height),
         zi_wT = as.numeric(zi_wT), zi_uw = as.numeric(zi_uw)) %>%
  pivot_longer(cols = c(uw, wT), names_to = "name", values_to = "value") %>%
  mutate(zi = ifelse(name == "uw", zi_uw, zi_wT)) %>%
  select(-c(zi_wT, zi_uw))

g1=
  ggplot(data_zi, aes(x=height,y=value, color=name, fill=name))+
  geom_point() + geom_line()+
  geom_vline(xintercept = 35,lty=2) +
  geom_hline(yintercept = 0.05,lty=2, color="grey") +
  geom_vline(aes(xintercept = zi, color=name),lty=2) +
  facet_wrap(wind_pred ~ stability, ncol = 6) +
  scale_x_continuous(breaks=c(unique(data_zi$height))) +
  coord_flip() +
  theme(strip.text = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text=element_text(size=12),
        legend.text = element_text(size = 12),
        axis.title=element_text(size=12))
g1

# correcting estimates
data_zi <- within(data_zi, zi[stability == "d.WS-2" & wind_pred == "b.90-180" & name == "uw"] <- 127)
data_zi <- within(data_zi, zi[stability == "f.VS" & wind_pred == "b.90-180" & name == "wT"] <- 100)
rm(g1, df_raw)


################### Apply the Lenschow model (1-z/h)^r ###########

sub_Lenschow = data_zi %>%
  mutate(height = as.numeric(height),
         zi = as.numeric(zi))

# Get the variables: h e r
sub_Lenschow_fit = sub_Lenschow %>%
  group_by(wind_pred,stability,name) %>%
  nest() %>%
  mutate(fit_Lenschow = map(data, ~ fitflux_LM(.x$value,.x$height,h_canopy,.x$zi[1],name)),
         fit_Lenschow_teo = map(data, ~ fitflux_Lenschow(.x$value,.x$height,h_canopy,.x$zi[1],name))) %>%
  unnest(cols = c(data,fit_Lenschow, fit_Lenschow_teo))

zfit = seq(35,316,1)

# Using the model
sub_Lenschow_plot = sub_Lenschow_fit %>%
  filter(height == 35) %>%
  dplyr::select(wind_pred,stability,name, zi, h, r, h_c, h_teo, r_teo) %>%
  group_by(wind_pred,stability,name) %>%
  nest() %>%
  mutate(height = map(data, ~ zfit),
         fit_Lenschow = map(data, ~ (1-(zfit-.x$h_c)/.x$h)^.x$r),
         fit_Lenschow_teo = map(data, ~ (1-(zfit-.x$h_c)/.x$h_teo)^.x$r_teo)) %>%
  unnest(cols = c(data, height, fit_Lenschow, fit_Lenschow_teo)) %>%
  mutate(z_on_h = (height - h_c)/h,
         z_on_hteo = (height - h_c)/h_teo)

# Corrigindo o fluxo
sub_Lenschow_correction = sub_Lenschow %>%
  group_by(wind_pred,stability,name) %>%
  nest() %>%
  mutate(correct_flux = map(data, ~ correct_norm(.x$value,.x$height))) %>%
  unnest(cols = c(data,correct_flux)) %>%
  dplyr::select(!c(value)) %>%
  rename(value = correct_flux) %>%
  mutate(z_on_h = height / zi,
         value2 = ifelse(value < 0, 10^(-3), value))



################### Figure 3: Lenschow ##########

# Dataframe com os valores de zi e zi by lenschow (at 5%)

data_zi = sub_Lenschow_plot %>% 
  group_by(wind_pred,stability,name) %>% 
  summarise(h_c = unique(h_c),
            zi = unique(zi),
            zi1 = getZiLenschow(fit_Lenschow, height),
            zi_teo = getZiLenschow(fit_Lenschow_teo, height)) 

names = c("wT"=expression(paste(bar("w'T'"))), 
          "uw"=expression(paste(bar("u'w'"))) ) 
colr = c("wT"="#E41A1C", "uw"="#2166AC")
winddir = c("a.00-70" = "Lower roughness", "b.90-180" = "Higher roughness")

gg_fig3 = 
  ggplot(sub_Lenschow_correction, 
         aes(y=value, x=height, color = name,fill = name)) +
  geom_point(size=2.25, shape=16) + 
  geom_point(aes(y=value2, x=height, color = name,fill = name),size=2, shape=10) +
  geom_smooth(se=F, linewidth = 1) +  
  geom_line(data = sub_Lenschow_plot, aes(x = height, y = fit_Lenschow_teo, color = name),lty=4) + 
  geom_hline(yintercept = 0.06, lty=2, color="black") +
  scale_y_log10(labels=trans_format('log10',math_format(10^.x)),
                limits = c(10^-3,3)) +
  scale_x_log10(labels=trans_format('log10',math_format(10^.x)),
                limits = c(3*10^1,4*10^2), 
                breaks = as.numeric(1 %o% 10 ^ (1:3)),
                minor_breaks = as.numeric(1:10 %o% 10 ^ (1:3))) +
  annotation_logticks(sides="tblr") +
  scale_color_manual(name = "Profile", values = colr, labels = names) +
  scale_fill_manual(name = "Profile", values = colr, labels = names) +
  coord_flip() + theme_bw() +
  theme(strip.text = element_text(size=14, face="bold"),
        strip.background = element_rect(colour="white", fill="white"),
        axis.text = element_text(size = 14, color = "black"), 
        plot.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        legend.background = element_rect( colour = "grey"),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14)) +
  labs(x = "Height (m)", y= "Normalized flux") +
  facet_wrap( wind_pred ~ stability, ncol = 6,
              labeller = labeller(wind_pred = as_labeller(winddir))) 
gg_fig3

ggsave(plot = gg_fig3,
       filename = paste("GRL/fig3.png",sep=""),
       width = 14, height = 9, units = "in", scale = 1, dpi = 300)


rm(names, colr,winddir,gg_fig3)



########################### Figure 4 ################

sub_Lenschow_table = data_zi %>% 
  dplyr::select(wind_pred,stability,name, zi, zi_teo, h_c)  %>%
  pivot_longer(!c(wind_pred,stability, name, h_c),
               names_to = "model")

winddir = c("a.00-70" = "a. Lower roughness", "b.90-180" = "b. Higher roughness")
colr = c("wT"="#E41A1C", "uw"="#2166AC")
names = c("wT"=expression(paste(bar("w'T'"))), 
          "uw"=expression(paste(bar("u'w'")))   )
names2 = c("zi"="Obs. data", "zi_teo"="Model")
eixox = c("a.NN" = "NN", "b.NN-WS" = "NN-WS", "c.WS-1" = "WS-1", 
          "d.WS-2" = "WS-2", "e.WS-VS" = "WS-VS", "f.VS" = "VS")
shp = c(zi = 25, zi_teo = 19)

gg_fig4 = 
  ggplot(sub_Lenschow_table,
         aes(x=stability,y = value, color = name, shape = model)) +
  geom_point(aes(size = model), stroke = 1.) +
  coord_cartesian(ylim = c(35,316)) +
  scale_size_manual(values = c(3,3)) + 
  scale_shape_manual(name = expression(h[N] ~ "estimated by"), values = shp, labels = names2) +
  scale_color_manual(name = "Profile", values = colr, labels = names) +
  scale_x_discrete(labels = eixox) + 
  scale_y_continuous(breaks = c(35,100,200,300)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16, color = "black"), 
        axis.title = element_text(size = 16), 
        strip.text = element_text(face="bold",size = 16),
        legend.text = element_text(size = 16), legend.title = element_text(size = 16),
        legend.background = element_rect( colour = "grey"),
        strip.background = element_rect(colour="white", fill="white"))+
  facet_grid( .~wind_pred , scale = "fixed",
              labeller = labeller(wind_pred = as_labeller(winddir)))+
  guides(size = "none")  +
  labs(x="Stability classes", y = "Height (m)")  
gg_fig4

ggsave(plot = gg_fig4,
       filename = paste("GRL/fig4.png",sep=""),
       width = 12, height = 6, units = "in", scale = 1, dpi=300)


rm(sub_Lenschow_table, winddir,colr,names,names2,eixox,shp,gg_fig4)





