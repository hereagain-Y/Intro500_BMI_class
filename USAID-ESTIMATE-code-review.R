# USAID ESTIMATE
library(biostat3)
library(tidyverse)
library(tsModel)
library(sandwich)
library(readxl)
library(stargazer)
library(dplyr)
library(tsModel)
library(sandwich)
library(lubridate)
library(ggthemes)
library(zoo)
library(readxl)
library(stargazer)
library(dplyr)
library(MASS)


# load data 
load("/Users/hereagain/JHU/RA_JHU/COVID19ITS/Updated_data/anova_msm.RData")
load("/Users/hereagain/JHU/RA_JHU/COVID19ITS/Updated_data/wrhi_fsw.RData")

# Newey adjustment to regression 
Newey_result <- function(model1){
  est <- exp(c(
    coef(model1)["High"], 
    coef(model1)["Low"], 
    coef(model1)["time"], 
    coef(model1)["High:inter1"] + coef(model1)["time"],
    coef(model1)["Low:inter2"] + coef(model1)["time"]
  ))
  # For adjusted 
  # for se 
  se1 <- sqrt(diag(NeweyWest(model1, prewhite = F, lag = 3)))["High"]
  se2 <- sqrt(diag(NeweyWest(model1, prewhite = F, lag = 3)))["Low"]
  se3 <- sqrt(diag(NeweyWest(model1, prewhite = F, lag = 3)))["time"]
  
  # lower 95%CI
  # ---------------------------------------------------------------------
  # CODE REVIEW: 
  # I generally break long lines like this into intermediate variables so 
  # it's more readable. I've assigned the vars lincom_lb_high_inter1 and 
  # lincom_lb_low_inter2 but these should be renamed to be informative.
  # ---------------------------------------------------------------------
  lb <- round(est[c(1:3)] * exp(-1.96 * c(se1, se2, se3)),2)
  lincom_lb_high_inter1 <- lincom(model = model1, c("time + High:inter1"), eform = T)[2]
  lincom_lb_low_inter2 <- lincom(model = model1, c("time + Low:inter2"), eform = T)[2]
  trend_lb <- as.numeric(c(lincom_lb_high_inter1, lincom_lb_low_inter2))
  lb<- round(c(lb, trend_lb), 2)
  
  # upper 95% CI
  # ---------------------------------------------------------------------
  # CODE REVIEW: 
  # Same feedback here, I'm not familiar with the lincom function 
  # but this is more readable for me with intermediate variables.
  # ---------------------------------------------------------------------
  ub <- est [c(1:3)]* exp(1.96 * c(se1, se2, se3))
  lincom_ub_high_inter1 <- lincom(model = model1, c("time + High:inter1"), eform = T)[3]
  lincom_ub_low_inter2 <- lincom(model = model1, c("time + Low:inter2"), eform = T)[3]
  trend_ub <- as.numeric(c(lincom_ub_high_inter1, lincom_ub_low_inter2))
  ub<- round(c(ub, trend_ub), 2)
  est2 <- round(est, 2)
  
  # ---------------------------------------------------------------------
  # CODE REVIEW:
  # I think the code here looks great but it might benefit from some
  # comments explaining what's being done.
  #
  # I put the following code as an example of a way I personally find more
  # readable. I like when I have an example string that variables get 
  # inserted into versus string chunks getting appended (like f strings in python).
  # I think you can also just do this with sprintf depending on your inputs:
  # i.e. irr1 <- sprintf("%s (%s, %s)", est2[1], lb[1], ub[1])
  #
  # This is just personal preference though and your approach is much more
  # concise!
  # ---------------------------------------------------------------------
  # EXAMPLE CODE:
  # pattern_str <- "%s (%s, %s)"
  # 
  # pattern_vals_irr1 <- c(est2[1], lb[1], ub[1])
  # pattern_vals_irr2 <- c(est2[2], lb[2], ub[2])
  # pattern_vals_beta_1 <- c(est2[3], lb[3], ub[3])
  # pattern_vals_trend1 <- c(est2[4], lb[4], ub[4])
  # pattern_vals_trend2 <- c(est2[5], lb[5], ub[5])
  # 
  # irr1 <- do.call(sprintf, as.list(c(pattern_str, pattern_vals_irr1)))
  # irr2 <- do.call(sprintf, as.list(c(pattern_str, pattern_vals_irr2)))
  # beta_1 <- do.call(sprintf, as.list(c(pattern_str, pattern_vals_beta_1)))
  # trend1 <- do.call(sprintf, as.list(c(pattern_str, pattern_vals_trend1)))
  # trend2 <- do.call(sprintf, as.list(c(pattern_str, pattern_vals_trend2)))
  # ---------------------------------------------------------------------
  
  irr1 = paste(est2[1]," (",lb[1],",",ub[1],")")
  irr2 = paste(est2[2]," (",lb[2],",",ub[2],")")
  beta_1 =  paste(est2[3]," (",lb[3],",",ub[3],")")
  
  trend1 = paste(est2[4]," (",lb[4],",",ub[4],")")
  trend2 = paste(est2[5]," (",lb[5],",",ub[5],")")
  
  result = as.data.frame(cbind(irr1, irr2, beta_1, trend1, trend2))
  names(result) <- c(
    "IRR1", 
    "IRR2", 
    "Pre-lock down Trend", 
    "Post lock down1 Trend", 
    "Post lock down2 Trend"
  )
  return(result)
}

# ---------------------------------------------------------------------
# CODE REVIEW:
# Add a comment here explaining this code, and maybe space it out to separate
# it from the function definition statements.
# ---------------------------------------------------------------------
exp(coef(model1_init)["High:inter1"] + coef(model1_init)["time"])

combine_nwresults = function(mod1, mod2, mod4, mod5){
  hiv = Newey_result(mod1)
  pos = Newey_result(mod2)
  prep = Newey_result(mod4)
  tx = Newey_result(mod5)
  result <- rbind(hiv, pos, prep, tx)
  row.names(result)<- c("HIV test",'POS tests',"Prep","ART")
  return(result)
}

# ---------------------------------------------------------------------
# CODE REVIEW:
# I think the rest of the code looks great but would definitely benefit from 
# some comments, without any background info it's difficult to interpret. 
# I reformatted it a little to break up the long lines but that's just 
# personal preference.
# ---------------------------------------------------------------------
model1_init  = glm.nb(
  HTS_TST_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_msm
)
nb1 =  glm.nb(
  HTS_TST_POS_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_msm
)
model3_init  = glm.nb(
  PREP_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_msm
)
tx1 = glm.nb(
  TX_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_msm
)

model_seas1 =  glm.nb(
  HTS_TST_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 37), 
  sum_msm
)
nb_seas1 =  glm.nb(
  HTS_TST_POS_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 37), 
  sum_msm
)
model_seas3 =  glm.nb(
  PREP_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2+harmonic(time, 2, 36), 
  sum_msm
)
tx_seas2 =  glm.nb(
  TX_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 37) + offset(log(HTS_TST_POS_TOTAL)), 
  sum_msm
)

fsw_test1  = glm.nb(
  HTS_TST_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_fsw
)
fsw_pos1  = glm(
  HTS_TST_POS_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  family = "poisson", 
  sum_fsw
)
fsw_prep2  = glm.nb(
  PREP_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_fsw
)
fsw_art2  = glm.nb(
  TX_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2, 
  sum_fsw
)

fsw_seas1 = glm.nb(
  HTS_TST_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 51), 
  sum_fsw
)
fsw_sea2 = glm(
  HTS_TST_POS_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 53), 
  family = "poisson",  
  sum_fsw
)
fsw_nb_sea3 = glm.nb(
  PREP_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 53), 
  sum_fsw
)
fsw_nb_sea4 = glm.nb(
  TX_NEW_TOTAL~time + High + Low + High:inter1 + Low:inter2 + harmonic(time, 2, 52), 
  sum_fsw
)

msm_re1 = combine_nwresults(model1_init, nb1, model3_init, tx1)
msm_seas = combine_nwresults(model_seas1, nb_seas1, model_seas3, tx_seas2)

fsw_re1 = combine_nwresults(fsw_test1, fsw_pos1, fsw_prep2, fsw_art2)
fsw_seas = combine_nwresults(fsw_seas1, fsw_sea2, fsw_nb_sea3, fsw_nb_sea4)

unadjusted_result = rbind(msm_re1, fsw_re1) 

adjusted_result = rbind(msm_seas, fsw_seas)

# export to a result table 
library(xlsx)

write.xlsx(
  unadjusted_result, 
  "/Users/hereagain/JHU/RA_JHU/COVID19ITS/USAID_NWadjuested_Summary_330_v3.xlsx", 
  sheetName = "Unadjusted", 
  col.names = TRUE, 
  row.names = TRUE, 
  append = FALSE, 
  showNA = FALSE
)
write.xlsx(
  adjusted_result, 
  "/Users/hereagain/JHU/RA_JHU/COVID19ITS/USAID_NWadjuested_Summary_330_v3.xlsx", 
  sheetName = "Adjusted", 
  col.names = TRUE, 
  row.names = TRUE, 
  showNA = FALSE, 
  append=TRUE
)
