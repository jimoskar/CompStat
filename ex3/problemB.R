setwd("~/Github/CompStat/ex3") # For jim

bilirubin <- read.table("files/bilirubin.txt",header=T)
head(bilirubin)


## 1. ----

# Create box plot of logarithm of measurement
ggplot(bilirubin, aes(x = pers, y = log(meas))) + 
  geom_boxplot() + xlab("pers") + ylab("log(meas)") +
  ggtitle("Log of Measured Bilirubin for Each Individual") + theme_minimal()
ggsave("figures/boxplot.pdf", height = 4, width = 6)


mod <- lm(log(meas) ~ pers, data = bilirubin) 
s <- summary(mod)
Fval <- s$fstatistic[1]

## 2. ----

permTest <- function(){
  df.perm <- data.frame(bilirubin)
  df.perm$pers <- sample(bilirubin$pers, size = nrow(bilirubin), 
                         replace = FALSE)
  mod <- lm(log(meas) ~ pers, data = df.perm)
  s <- summary(mod)
  return(s$fstatistic[1])
}

n <- 999 # Number of samples
F.samples <- rep(NA, n)
set.seed(4300)
for(i in 1:n){
  F.samples[i] <- permTest()
}

sum(F.samples >= Fval)/n # p-value
