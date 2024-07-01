# ==============================================================================
# load resources
# ==============================================================================

library(xts)
library(quantmod)
library(fredr)
fredr_set_key("32ff0f847aab35def1d3d5f7216e4064")

# ==============================================================================
# load old CCM (2016 JBES) data
# ==============================================================================

ccm_2016_JBES_14var_1964Q1_2013Q4 <- read.csv("~/DynamicCopula/datasets/ccm_2016_JBES_14var_1964Q1_2013Q4.csv", header = FALSE)

var_names <- c("gdpgrowth", "pcegrowth", "bfigrowth", "resinvgrowth", "ipgrowth", "capu", 
               "employ", "hours", "ur", "gdppigrowth", "pcepigrowth", "ffr", "spread", "realSPgrowth")

fred_handles <- c("GDPC1", "PCE", "PNFI", "PRFI", "INDPRO", "CUMFNS", "PAYEMS",
                  "AWHI", "UNRATE", "GDPCTPI", "PCEPI", "DFF")

n <- ncol(ccm_2016_JBES_14var_1964Q1_2013Q4)
T <- nrow(ccm_2016_JBES_14var_1964Q1_2013Q4)
start_date <- as.Date("1964-01-01")
end_date <- as.Date("2013-12-31")
dates = seq(1964.00, 2013.75, by = 0.25)

# ==============================================================================
# Can I perfectly replicate the CCM (2016 JBES) data (modulo subsequent revisions)
# ==============================================================================

level <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)

par(mfrow = c(n / 2, 2), mar = c(2, 2, 2, 2))

for(i in 1:n){
  
  if(var_names[i] == "spread"){
    bondyield <- fredr(
      series_id = "DGS10",
      observation_start = start_date,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    billyield <- fredr(
      series_id = "DTB3",
      observation_start = start_date,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    raw_data = bondyield - billyield
    
  }else if(var_names[i] == "realSPgrowth"){
    getSymbols("^GSPC", src = "yahoo", from = start_date - 365 / 4, to = end_date)
    
    # Aggregate the data to a quarterly frequency by averaging
    GSPC_quarterly <- apply.quarterly(GSPC, colMeans, na.rm = TRUE)
    
    #sp <- as.vector(GSPC_quarterly$GSPC.Close[145:345]) # AUTOMATE THIS!
    sp <- as.vector(GSPC_quarterly$GSPC.Close)
    
    pcepi <- fredr(
      series_id = "PCEPI",
      observation_start = start_date - 365 / 4,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    raw_data <- sp / pcepi
  }else{
    raw_data <- fredr(series_id = fred_handles[i],
                      observation_start = start_date - !level[i]*365/4,
                      observation_end = end_date,
                      frequency = "q", 
                      aggregation_method = "avg")$value
  }
  
  if(level[i] == TRUE){
    freddat = raw_data
  }else{
    freddat = 400 * diff(log(raw_data))
  }
  
  plot(dates, ccm_2016_JBES_14var_1964Q1_2013Q4[, i], type = 'l', main = var_names[i])
  
  if(var_names[i] == "hours"){
    lines(dates[2:T], freddat, col = "red")
  }else{
    lines(dates, freddat, col = "red")
  }
  
}

# ==============================================================================
# Can I construct an updated series?
# ==============================================================================

start_date <- as.Date("1965-01-01")
end_date <- as.Date("2019-12-31")
dates <- seq(1965.00, 2019.75, by = 0.25)

T <- length(dates)

Y <- matrix(0, T, n)

for(i in 1:n){
  
  if(var_names[i] == "spread"){
    bondyield <- fredr(
      series_id = "DGS10",
      observation_start = start_date,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    billyield <- fredr(
      series_id = "DTB3",
      observation_start = start_date,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    raw_data = bondyield - billyield
    
  }else if(var_names[i] == "realSPgrowth"){
    getSymbols("^GSPC", src = "yahoo", from = start_date - 365 / 4, to = end_date)
    
    # Aggregate the data to a quarterly frequency by averaging
    GSPC_quarterly <- apply.quarterly(GSPC, colMeans, na.rm = TRUE)
    
    #sp <- as.vector(GSPC_quarterly$GSPC.Close[145:345]) # AUTOMATE THIS!
    sp <- as.vector(GSPC_quarterly$GSPC.Close)
    
    pcepi <- fredr(
      series_id = "PCEPI",
      observation_start = start_date - 365 / 4,
      observation_end = end_date,
      frequency = "q", 
      aggregation_method = "avg"
    )$value
    
    raw_data <- sp / pcepi
  }else{
    raw_data <- fredr(series_id = fred_handles[i],
                      observation_start = start_date - !level[i]*365/4,
                      observation_end = end_date,
                      frequency = "q", 
                      aggregation_method = "avg")$value
  }
  
  if(level[i] == TRUE){
    Y[, i] = raw_data
  }else{
    Y[, i] = 400 * diff(log(raw_data))
  }
  #plot(dates, Y[, i], type = 'l', main = var_names[i])
  
}

# ==============================================================================
# plot it!
# ==============================================================================

my_labels <- c("Real GDP", "PCE", "BFI", "Residential invesment",
               "Industrial production", "Capacity utilization", "Payroll employment",
               "Hours", "Unemployment rate", "GDP PI inflation", "PCE inflation",
               "Federal funds rate", "Yield spread", "Real S&P 500")

if(write_image){
  png(paste("_images/ccm_2016_JBES_14var_1965Q1_2019Q4.png", sep = ""), 
      width = 6.5, height = 8, units = "in", res = 650)
}

write_image <- TRUE
histogram <- FALSE
add_normal <- FALSE

par(mfrow = c(n / 2, 4))

for(i in 1:n){
  y = Y[, i]
  
  par(mar = c(2, 2, 2, 0.25))
  plot(dates, y, type = "l", main = my_labels[i])
  
  par(mar = c(2, 0.25, 2, 0.25))
  if(histogram == TRUE){
    hist(y, breaks = 25, freq = FALSE, col = "blue", border = "white", main = "", yaxt = "n")
  }else{
    x_grid <- seq(min(y)-5, max(y)+5, length.out = 750)
    y_vals <- dnorm(x_grid, mean = mean(y), sd = sd(y))
    f = density(y)
    
    plot(f, yaxt = "n", bty = "n", main = "", col = "blue", ylim = c(0, max(max(f$y), max(y_vals))))
    if(add_normal == TRUE){
      lines(x_grid, y_vals, col = "red")
    }
  }
}

if(write_image){
  dev.off()
}








