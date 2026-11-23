full_nsw_crime_dataset <- read.csv("~/DynamicCopula/datasets/full_nsw_crime_dataset.csv")

T <- nrow(full_nsw_crime_dataset)
n <- ncol(full_nsw_crime_dataset) - 1

par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))

for(i in 1:n){
  y = full_nsw_crime_dataset[, i + 1]
  minY = min(y)
  maxY = max(y)
  plot(1:T, y, type = "l", main = colnames(full_nsw_crime_dataset)[i + 1])
  probs = table(factor(y, levels = minY:maxY)) / T
  if(length(unique(y)) <= 100){
    barplot(probs, names.arg = minY:maxY)
  } else {
    hist(y, breaks = "Scott", freq = FALSE)
  }
}