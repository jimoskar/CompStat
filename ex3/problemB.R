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
  df.perm$pers <- sample(bilirubin$pers, replace = FALSE)
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

# Create histogram
F.right <- F.samples[which(F.samples >= Fval)]
F.left <- setdiff(F.samples, F.right)
F.df <- data.frame(F.vals = c(F.left, F.right), 
                   ind = as.factor(c(rep( "< Fval", length(F.left)),rep( "> Fval", length(F.right)))))
ggplot(F.df) + 
  geom_histogram(aes(x = F.vals, colour = ind), binwidth = 0.1, fill = "white") +
  xlab("F-values") + ylab("Count") +
  theme_minimal() +  theme(legend.title=element_blank())
ggsave("figures/F-vals.pdf", width = 7, height = 4)
