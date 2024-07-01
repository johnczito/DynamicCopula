Y <- read.csv("datasets/sleep_data.csv")

T <- nrow(Y)
n <- ncol(Y)
labels = c("Heart rate", "EEG sleep state", "Body temperature")
states = c("quiet", "?", "active", "awake")

par(mfcol = c(2, n), mar = c(2, 2, 2, 2))

for(i in 1:n){

  plot(1:T, Y[, i], type = "l", main = labels[i], cex.main = 2,
       yaxt = "n")
  if(i == 2){
    axis(2, at = 1:4, labels = states, las = 1)
  }else{
    axis(2)
  }
  
  if(i == 1){
    #plot(density(Y[, i]), main = "")
    hist(Y[, i], breaks = "Scott", main = "", col = "lightblue")
  }else if(i == 2){
    barplot(table(Y[, i]) / T, names.arg = states, col = "lightblue")
  }else{
    barplot(table(Y[, i]) / T, names.arg = sort(unique(Y[, i])), col = "lightblue")
  }
}