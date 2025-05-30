rm(list = ls())
library(foreign)
data.ch3 <- read.dta(file = "D:/Desktop/Stat Comput-Karlis/racd03data.dta")
attach(data.ch3)   
head(data.ch3)
X = data.ch3$DVISITS
Y = data.ch3$MEDICINE
# Table 1
t1 = table(X, Y)
t1

