setwd("~/Github/CompStat/ex3") # For jim

bilirubin <- read.table("data/bilirubin.txt",header=T)
head(bilirubin)


## 1. ----

# Create box plot of logarithm of measurement
ggplot(bilirubin, aes(x = pers, y = log(meas))) + 
   geom_boxplot() + xlab("Person") + ylab("log(Measurement)") + theme_minimal()
ggsave("figures/boxplot.pdf", height = 5, width = 5)


mod <- lm(log(meas) ~ -1 +pers, data = bilirubin) 
s <- summary(mod)
F.val <- s$fstatistic[1]

## 2. ----

permTest <- function(){
  perm <- sample(c("p1", "p2", "p3"), size = nrow(bilirubin), replace = TRUE)
  df.perm <- data.frame(meas = bilirubin$meas, pers = perm)
  mod <- lm(log(meas) ~ -1 +pers, data = df.perm)
  s <- summary(mod)
  return(s$fstatistic[1])
}

permTest()

F.vals <- rep(NA, 1000)
F.vals[1] = F.val
for(i in 2:1000){
  F.vals[i] <- permTest()
}

sum(F.vals >= F.val)/length(F.vals)
