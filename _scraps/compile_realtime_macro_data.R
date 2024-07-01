var_names <- c("gdpgrowth", "pcegrowth", "bfigrowth", "resinvgrowth", "ipgrowth", "capu", 
               "employ", "hours", "ur", "gdppigrowth", "pcepigrowth", "ffr", "spread", "realSPgrowth")

unrevised <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_ur_ffr_spread_1964Q1_2014Q1.csv", header = FALSE)
bfigrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_bfigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
capu <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_capu_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
employ <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_employ_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdpgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_gdpgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
gdppigrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_gdppigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
hours <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_hours_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
ipgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_ipgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcegrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_pcegrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
pcepigrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_pcepigrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
realSPgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_realSPgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)
resinvgrowth <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_resinvgrowth_vintages_1985Q1_2014Q2_backto_1964Q1.csv", header = FALSE)

dates = seq(1964.00, 2014.00, by = 0.25)
vintages <- seq(1985.00, 2014.25, by = 0.25)

n <- 14
T <- nrow(unrevised)
V <- length(vintages)

var_names <- c("gdpgrowth", "pcegrowth", "bfigrowth", "resinvgrowth", "ipgrowth", "capu", 
               "employ", "hours", "ur", "gdppigrowth", "pcepigrowth", "ffr", "spread", "realSPgrowth")

real_time_macro_data <- array(0, c(T, n, V))

for(v in 1:V){
  vintage <- vintages[v]
  end_period <- vintage - 0.25
  T <- which(dates == end_period)
  
  real_time_macro_data[1:T, 01, v] <- gdpgrowth[1:T, v]
  real_time_macro_data[1:T, 02, v] <- pcegrowth[1:T, v]
  real_time_macro_data[1:T, 03, v] <- bfigrowth[1:T, v]
  real_time_macro_data[1:T, 04, v] <- resinvgrowth[1:T, v]
  real_time_macro_data[1:T, 05, v] <- ipgrowth[1:T, v]
  real_time_macro_data[1:T, 06, v] <- capu[1:T, v]
  real_time_macro_data[1:T, 07, v] <- employ[1:T, v]
  real_time_macro_data[1:T, 08, v] <- hours[1:T, v]
  real_time_macro_data[1:T, 09, v] <- unrevised[1:T, 1]
  real_time_macro_data[1:T, 10, v] <- gdppigrowth[1:T, v]
  real_time_macro_data[1:T, 11, v] <- pcepigrowth[1:T, v]
  real_time_macro_data[1:T, 12, v] <- unrevised[1:T, 2]
  real_time_macro_data[1:T, 13, v] <- unrevised[1:T, 3]
  real_time_macro_data[1:T, 14, v] <- realSPgrowth[1:T, v]
}