# Plots



library(dplyr)
# library(readstata13)
# library(xlsx)
source("H:/My R functions/MySummary.R")
source("H:/My R functions/lintestOR.R")
source("H:/My R functions/tidyoutput.R")
source("H:/My R functions/niceN.R")
source("H:/My R functions/niceP.R")
# library(janitor)
# library(officer)
# library(flextable)
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(ggplot2)
library(survsup)
# library(epitools)
library(psych)
library(lme4)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(car)
library(extrafont)
loadfonts()
fonts()
library(forcats)


dat <- readRDS("Z:/Group_work/Alex/Encrypted/Alex/Adult Asthma/SCC 2019-2020/Data/tidyData/AA_SCC_2019-20_clean_data_2020-06-23.RDS")


hist(dat$heart_rate)
hist(dat$resp_rate)


# For arrival to beta agonists:

# # KM curve
# # We don't want Scotland anymore
# 
# datKMPEF <- filter(dat, !is.na(arrivaltoPEF)) %>% filter(country != "Scotland")
# datKMPEF$seen <- 1
# 
# # Optionally we can replicate this to give all, but it basically just follows the English curve
# # sccKMall <- sccKM
# # sccKMall$country <- "All"
# # Now we bind it together
# # sccKM <- rbind(sccKM, sccKMall)
# 
# head(sort(datKMPEF$arrivaltoPEF))
# 
# datKMPEF %>% filter(country != "Scotland") %>% filter(arrivaltoPEF<(-0.0001)) %>% select(arrivaltoPEF) %>% nrow()
# 
# # If they arrive earlier then just change it to zero for the sake of the kaplan-meier 
# datKMPEF$arrivaltoPEF[datKMPEF$arrivaltoPEF<(-0.000001)] <- 0
# 
# 
# head(sort(datKMPEF$arrivaltoPEF), 10)
# 
# # 48 hours
# survfit(Surv(datKMPEF$arrivaltoPEF*24, datKMPEF$seen) ~ datKMPEF$country, data = datKMPEF) %>% 
#   plot_survfit(ci = TRUE, legend.title = "Country", xmax = 48, xbreaks = seq(0, 48, 6)) + 
#   labs(x = "Time (hours)", y = "Percentage of patients who have received PEF (%)")
# 
# 
# 
# 
# 





# 3.4.3 KM curve
# remove Scotland

dat$arrival_to_b2a_minutes %>% head()

datKMBA <- filter(dat, !is.na(arrival_to_b2a_minutes)) # %>% filter(country != "Scotland")
datKMBA$seen <- 1

# Optionally we can replicate this to give all, but it basically just follows the English curve
# sccKMall <- sccKM
# sccKMall$country <- "All"
# Now we bind it together
# sccKM <- rbind(sccKM, sccKMall)


# in minutes
# survfit(Surv(datKMBA$arrivaltob2agonists*24*60, datKMBA$seen) ~ datKMBA$country, data = datKMBA) %>% 
#   plot_survfit(ci = TRUE, legend.title = "Country", xmax = 360, xbreaks = seq(0, 360, 30)) + 
#   labs(x = "Time (minutes)", y = "Percentage of patients who have received PEF within 48 hours (%)")

# in hours
survfit(Surv(datKMBA$arrival_to_b2a_minutes/60, datKMBA$seen) ~ datKMBA$country, data = datKMBA) %>% 
  plot_survfit(ci = TRUE, legend.title = "Country", xmax = 12, xbreaks = seq(0, 12, 1)) + 
  labs(x = "Time (hours)", 
       y = expression(paste("Percentage of patients who have received ", beta[2]," agonists (%)")))



# And steroids:

dat$arrival_to_steroids_hours %>% head()

datKMS <- filter(dat, !is.na(arrival_to_steroids_hours)) # %>% filter(country != "Scotland")
datKMS$seen <- 1

# Optionally we can replicate this to give all, but it basically just follows the English curve
# sccKMall <- sccKM
# sccKMall$country <- "All"
# Now we bind it together
# sccKM <- rbind(sccKM, sccKMall)


# in minutes
# survfit(Surv(datKMS$arrivaltob2agonists*24*60, datKMS$seen) ~ datKMS$country, data = datKMS) %>% 
#   plot_survfit(ci = TRUE, legend.title = "Country", xmax = 360, xbreaks = seq(0, 360, 30)) + 
#   labs(x = "Time (minutes)", y = "Percentage of patients who have received PEF within 48 hours (%)")

# in hours
survfit(Surv(datKMS$arrival_to_steroids_hours, datKMS$seen) ~ datKMS$country, data = datKMS) %>% 
  plot_survfit(ci = TRUE, legend.title = "Country", xmax = 48, xbreaks = seq(0, 48, 4)) + 
  labs(x = "Time (hours)", 
       y = "Percentage of patients who have received systemic steroids (%)")





# And PEF


dat$arrival_to_PEF_init_hours %>% head()

datKMPEF <- filter(dat, !is.na(arrival_to_PEF_init_hours)) # %>% filter(country != "Scotland")
datKMPEF$seen <- 1
datKMPEF$arrival_to_PEF_init_hours[datKMPEF$arrival_to_PEF_init_hours < 0] <- 0


# Optionally we can replicate this to give all, but it basically just follows the English curve
# sccKMall <- sccKM
# sccKMall$country <- "All"
# Now we bind it together
# sccKM <- rbind(sccKM, sccKMall)


# in minutes
# survfit(Surv(datKMPEF$arrivaltob2agonists*24*60, datKMPEF$seen) ~ datKMPEF$country, data = datKMPEF) %>% 
#   plot_survfit(ci = TRUE, legend.title = "Country", xmax = 360, xbreaks = seq(0, 360, 30)) + 
#   labs(x = "Time (minutes)", y = "Percentage of patients who have received PEF within 48 hours (%)")

# in hours
survfit(Surv(datKMPEF$arrival_to_PEF_init_hours, datKMPEF$seen) ~ datKMPEF$country, data = datKMPEF) %>% 
  plot_survfit(ci = TRUE, legend.title = "Country", xmax = 48, xbreaks = seq(0, 48, 4)) + 
  labs(x = "Time (hours)", 
       y = "Percentage of patients who have received an initial peak flow reading")



