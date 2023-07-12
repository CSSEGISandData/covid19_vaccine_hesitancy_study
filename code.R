library(mgcv)
library(ggplot2)
library(rstanarm)
library(tidyr)
library(broom.mixed)
library(dplyr)
library(gratia)
library(dotwhisker)
library(gridExtra)

setwd(getwd())
dat <- read.csv("dataset.csv")

# NOTE: due to the redistribution rules of Nielsen, the cable TV viewership data was excluded from the shared dataset.

##########################################################################
# GAM functions
##########################################################################
# weighted model
run_gam_func_wt <- function(dat_input) {
  gam_wt_output <- gam(VH_2021_12_15 ~  s(case_rate, k=3) + 
                         s(black_pct, k=3) + s(hispanic_pct, k=3) + 
                         s(higher_edu, k=3) +
                         s(med_hhincome, k=3) + 
                         s(median_age, k=3) + 
                         s(wo_insur_pct, k=3) +
                         s(vehicle_perhh, k=3) + 
                         s(std_fox_rtg, k=3) + s(std_cnn_rtg, k=3) + 
                         s(std_msn_rtg, k=3) + s(std_local_rtg, k=3) +
                         s(MMR_VR, k=3) + 
                         s(republican, k=3),
                       data = dat_input, method = "REML", weights=wt, 
                       select = TRUE, family = gaussian(link = log))
  
  return(gam_wt_output)
}

# unweighted model
run_gam_func_uw <- function(dat_input) {
  gam_uw_output <- gam(VH_2021_12_15 ~  s(case_rate, k=3) + 
                         s(black_pct, k=3) + s(hispanic_pct, k=3) + 
                         s(higher_edu, k=3) +
                         s(med_hhincome, k=3) + 
                         s(median_age, k=3) + 
                         s(wo_insur_pct, k=3) +
                         s(vehicle_perhh, k=3) + 
                         s(std_fox_rtg, k=3) + s(std_cnn_rtg, k=3) + 
                         s(std_msn_rtg, k=3) + s(std_local_rtg, k=3) +
                         s(MMR_VR, k=3) + 
                         s(republican, k=3),
                       data = dat_input, method = "REML",
                       select = TRUE, family = gaussian(link = log))
  
  return(gam_uw_output)
}


outliers <- function(x) {

  mean <- mean(x, na.rm=TRUE)
  std <- sd(x, na.rm=TRUE)
  upper_limit <- mean+(4*std)
  x > upper_limit
}

detect_outliers <- function(df, cols = names(df)) {
  for (col in cols) {
    df <- df[outliers(df[[col]]),]
  }
  df[!is.na(df$FIPS),]
}

# remove outliers
odf1 <- detect_outliers(dat, c('case_rate_12_15_21'))
odf2 <- detect_outliers(dat, c('std_cnn_rtg'))
odf3 <- detect_outliers(dat, c('std_fox_rtg'))
odf4 <- detect_outliers(dat, c('std_msn_rtg'))
odf5 <- detect_outliers(dat, c('std_local_rtg'))
outliers_dat <- Reduce(union, list(odf1,odf2,odf3,odf4,odf5))
dat <- dat[dat$FIPS %in% c(setdiff(dat$FIPS,outliers_dat$FIPS)),]

# normalize the weights
dat$wt <- log(dat$pop2020)/mean(log(dat$pop2020),na.rm=TRUE)

dat_1215 <- dat
names(dat_1215)[names(dat_1215) == "case_rate_12_15_21"] <- "case_rate"

##########################################################################
# GAMs - primary model

# NOTE: due to the redistribution rules of Nielsen, the cable TV viewership data was excluded from the shared dataset. This may affect the coding output.
##########################################################################

gam_w_dec <- gam(VH_2021_12_15 ~ s(case_rate, k=3) + 
                   s(black_pct, k=3) + s(hispanic_pct, k=3) + 
                   s(higher_edu, k=3) +
                   s(med_hhincome, k=3) + 
                   s(median_age, k=3) + 
                   s(wo_insur_pct, k=3) +
                   s(vehicle_perhh, k=3) + 
                   s(std_fox_rtg, k=3) + s(std_cnn_rtg, k=3) +
                   s(std_msn_rtg, k=3) + s(std_local_rtg, k=3) +
                   s(MMR_VR, k=3) + 
                   s(republican, k=3),
               data = dat_1215, method = "REML", weights=wt, 
               select = TRUE, family = gaussian(link = "log"))

draw(gam_w_dec)
summary(gam_w_dec)
concurvity(gam_w_dec)

##########################################################################
# GAMs - Sub-model incorporating Twitter misinformation rates 
##########################################################################

gam_w_dec <- gam(VH_2021_12_15 ~ s(case_rate, k=3) + 
                   s(black_pct, k=3) + s(hispanic_pct, k=3) + 
                   s(higher_edu, k=3) +
                   s(med_hhincome, k=3) + 
                   s(median_age, k=3) + 
                   s(wo_insur_pct, k=3) +
                   s(vehicle_perhh, k=3) + 
                   s(std_fox_rtg, k=3) + s(std_cnn_rtg, k=3) +
                   s(std_msn_rtg, k=3) + s(std_local_rtg, k=3) +
                   s(MMR_VR, k=3) + 
                   s(republican, k=3) +
                   s(twt501, k=3),
                 data = dat_1215, method = "REML", weights=wt, 
                 select = TRUE, family = gaussian(link = "log"))

draw(gam_w_dec)
summary(gam_w_dec)
concurvity(gam_w_dec)

##########################################################################
# GAMs - Land-use cluster-based sensitivity analysis 
##########################################################################

counties_dat_urban <- dat_1215[dat_1215$Area_Type == 'Urban',]
counties_dat_rural <- dat_1215[dat_1215$Area_Type == 'Rural',]
counties_dat_rural <- counties_dat_rural[!is.na(counties_dat_rural$FIPS),]
counties_dat_urban <- counties_dat_urban[!is.na(counties_dat_urban$FIPS),]

gam_urban <- run_gam_func_wt(counties_dat_urban)
gam_rural <- run_gam_func_wt(counties_dat_rural)
compare_smooths(gam_urban, gam_rural) %>% draw()
summary(gam_urban)
summary(gam_rural)


##########################################################################
# GAMs - clustered by population - no weights
##########################################################################

counties_dat_pop1 <- dat_1215[dat_1215$pop_cat == 'Q1',]
counties_dat_pop2 <- dat_1215[dat_1215$pop_cat == 'Q2',]
counties_dat_pop3 <- dat_1215[dat_1215$pop_cat == 'Q3',]
counties_dat_pop4 <- dat_1215[dat_1215$pop_cat == 'Q4',]

# GAMs with no weights
gam_pop1 <- run_gam_func_uw(counties_dat_pop1)
gam_pop2 <- run_gam_func_uw(counties_dat_pop2)
gam_pop3 <- run_gam_func_uw(counties_dat_pop3)
gam_pop4 <- run_gam_func_uw(counties_dat_pop4)
compare_smooths(gam_pop1, gam_pop2, gam_pop3, gam_pop4) %>% draw()


##########################################################################
# GAMs - visualize the output - primary model
##########################################################################

s_gam_w_dec <- smooth_estimates(gam_w_dec) %>% add_confint()

plot_1_gam <- function(x_s, x_var, x_lab) {
  
  Data_prob <- s_gam_w_dec %>% filter(smooth == x_s)
  
  ggplot(Data_prob) +
    geom_line(aes_string(x = x_var, y = "est")) +
    geom_ribbon(aes_string(x = x_var, ymax="upper_ci", ymin="lower_ci"), alpha=.2) +
    geom_rug(data=gam_w_dec$model, aes_string(x = x_var), sides = "b", length = grid::unit(0.03, "npc"))+
    labs(x = x_lab, y = "", title = "") +
    # ylim(-1, 0.7) +
    theme_bw()
}

p1 <- plot_1_gam("s(case_rate)","case_rate", 'COVID-19 case rate') + theme(legend.position = "none")
p2 <- plot_1_gam("s(black_pct)","black_pct", "Black (%)") + theme(legend.position = "none")
p3 <- plot_1_gam("s(hispanic_pct)","hispanic_pct", 'Hispanic (%)') + theme(legend.position = "none")
p4 <- plot_1_gam("s(higher_edu)","higher_edu", 'Postsecondary education (%)') + theme(legend.position = "none")
p5 <- plot_1_gam("s(med_hhincome)","med_hhincome", 'Median household income') + theme(legend.position = "none")
p6 <- plot_1_gam("s(median_age)","median_age", 'Median age') + theme(legend.position = "none")
p7 <- plot_1_gam("s(wo_insur_pct)","wo_insur_pct", 'Uninsured (%)') + theme(legend.position = "none")
p8 <- plot_1_gam("s(vehicle_perhh)","vehicle_perhh", 'Vehicles per household') + theme(legend.position = "none")
p9 <- plot_1_gam("s(std_fox_rtg)","std_fox_rtg", 'FNC') + theme(legend.position = "none")
p10 <- plot_1_gam("s(std_cnn_rtg)","std_cnn_rtg", 'CNN') + theme(legend.position = "none")
p11 <- plot_1_gam("s(std_msn_rtg)","std_msn_rtg", 'MSNBC') + theme(legend.position = "none")
p12 <- plot_1_gam("s(std_local_rtg)","std_local_rtg", 'Local news') + theme(legend.position = "none")
p13 <- plot_1_gam("s(MMR_VR)","MMR_VR", 'MMR coverage') + theme(legend.position = "none")
p14 <- plot_1_gam("s(republican)","republican", 'Republican (%)') + theme(legend.position = "none")

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,
             ncol = 4, widths = c(1, 1, 1, 1), left = "Partial Effect")

