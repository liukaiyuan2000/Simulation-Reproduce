install.packages("MASS")
install.packages("fda.usc")
install.packages("latex2exp")
install.packages("ggplot2")
install.packages("cowplot")
library(MASS)
library(fda.usc)
library(latex2exp)
library(ggplot2)
library(cowplot)
library(tibble)

n=100

# resultfir4=c(72.85,31.51,32.02,32.49,35.02,36.33,38.98,40.02,41.91,42.86)
# resultfir8=c(78.12,40.06,45.21,48.24,48.62,55.69,57.13,62.94,66.40,71.68)
# resultsec4=c(35.53,32.70,32.77,35.55,37.08,37.61,39.77,41.33,42.32,43.82)
# resultsec8=c(44.25,43.38,47.61,50.05,55.52,56.51,61.43,67.21,74.18,80.41)
result28=c(51.76,50.19,55.60,56.19,63.78,68.04,73.88,83.29,84.24,98.53)
result24=c(37.29,33.91,35.33,36.46,38.78,40.00,40.53,43.65,42.98,46.48)


result38=c(149.52,140.97,154.59,198.08,200.69,231.14,258.23,260.21,271.76,265.85)
result34=c(49.99,47.43,48.92,50.98,52.96,57.12,60.16,61.27,66.20,67.61)

result18=c(67.41,44.89,49.44,52.95,56.06,69.63,64.17,70.97,76.78,114.58)
result14=c(73.41,33.31,34.86,35.45,39.32,40.67,42.94,44.78,45.82,48.06)



par(mfrow=c(3,2),mar=c(5,6,4,2)+0.1)


barplot(result18/n,xlab = "m",ylab = TeX('mean value of \\mathbf{CV}(m)'),names.arg = c(1:10),main=TeX('Setting 1,$\\Sigma_1$'),cex.main=1.5,cex.lab=1.2,ylim=c(0,1.2))
barplot(result14/n,xlab = "m",ylab = TeX('mean value of \\mathbf{CV}(m)'),names.arg = c(1:10),main=TeX('Setting 1,$\\Sigma_2$'),cex.main=1.5,cex.lab=1.2,ylim=c(0,1.2))
barplot(result28/n,xlab = "m",ylab = TeX('mean value of \\mathbf{CV}(m)'),names.arg = c(1:10),main=TeX('Setting 2,$\\Sigma_1$'),cex.main=1.5,cex.lab=1.2,ylim=c(0,1))
barplot(result24/n,xlab = "m",ylab = TeX('mean value of \\mathbf{CV}(m)'),names.arg = c(1:10),main=TeX('Setting 2,$\\Sigma_2$'),cex.main=1.5,cex.lab=1.2,ylim=c(0,1))
barplot(result38/n,xlab = "m",ylab = TeX('mean value of \\mathbf{CV}(m)'),names.arg = c(1:10),main=TeX('Setting 3,$\\Sigma_1$'),cex.main=1.5,cex.lab=1.2,ylim=c(0,3))
barplot(result34/n,xlab = "m",ylab = TeX('mean value of \\mathbf{CV}(m)'),names.arg = c(1:10),main=TeX('Setting 3,$\\Sigma_2$'),cex.main=1.5,cex.lab=1.2,ylim=c(0,1))



err4eg1=read.csv("E:/first/eg1-44-cv.xls")[,-1]
err8eg1=read.csv("E:/first/eg1-88-cv.xls")[,-1]



df41 <- tibble(
  label = c(1:10),
  meanv = apply(err4eg1/100,2,mean),
  lower = meanv - apply(err4eg1/100,2,sd),
  upper = meanv + apply(err4eg1/100,2,sd),
  #group = c(rep("Group-1",7),rep("Group-2",7),rep("Group-3",8))
)


p1 <- ggplot(data = df41)+geom_point(aes(x = label, y = meanv),size=3)+
  geom_errorbar(aes(x = label,ymin = lower, ymax = upper),
                width = 0.6, # 控制上下两条短横线的长短
                size = 1 # 控制线条整体粗细
  )+theme(plot.title=element_text(hjust=0.5))+
  labs(x="m",y=TeX('mean value of \\mathbf{CV}(m)'),title =TeX('Setting 1,$\\Sigma_2$'))+
  xlim(0.5,10.5)+ylim(0,1.3)
p1


df81 <- tibble(
  label = c(1:10),
  meanv = apply(err8eg1/100,2,mean),
  lower = meanv - apply(err8eg1/100,2,sd),
  upper = meanv + apply(err8eg1/100,2,sd),
  #group = c(rep("Group-1",7),rep("Group-2",7),rep("Group-3",8))
)


p2 <- ggplot(data = df81)+geom_point(aes(x = label, y = meanv),size=3)+
  geom_errorbar(aes(x = label,ymin = lower, ymax = upper),
                width = 0.6, # 控制上下两条短横线的长短
                size = 1 # 控制线条整体粗细
  )+theme(plot.title=element_text(hjust=0.5))+
  labs(x="m",y=TeX('mean value of \\mathbf{CV}(m)'),title =TeX('Setting 1,$\\Sigma_1$'))+
  xlim(0.5,10.5)+ylim(-0.2,1.8)
p2

err4eg2=read.csv("E:/first/eg2-44-cv.xls")[,-1]
err8eg2=read.csv("E:/first/eg2-88-cv.xls")[,-1]


df42 <- tibble(
  label = c(1:10),
  meanv = apply(err4eg2/100,2,mean),
  lower = meanv - apply(err4eg2/100,2,sd),
  upper = meanv + apply(err4eg2/100,2,sd),
  #group = c(rep("Group-1",7),rep("Group-2",7),rep("Group-3",8))
)


p3 <- ggplot(data = df42)+geom_point(aes(x = label, y = meanv),size=3)+
  geom_errorbar(aes(x = label,ymin = lower, ymax = upper),
                width = 0.6, # 控制上下两条短横线的长短
                size = 1 # 控制线条整体粗细
  )+theme(plot.title=element_text(hjust=0.5))+
  labs(x="m",y=TeX('mean value of \\mathbf{CV}(m)'),title =TeX('Setting 2,$\\Sigma_2$'))+
  xlim(0.5,10.5)+ylim(0,1)
p3


df82 <- tibble(
  label = c(1:10),
  meanv = apply(err8eg2/100,2,mean),
  lower = meanv - apply(err8eg2/100,2,sd),
  upper = meanv + apply(err8eg2/100,2,sd),
  #group = c(rep("Group-1",7),rep("Group-2",7),rep("Group-3",8))
)



p4 <- ggplot(data = df82)+geom_point(aes(x = label, y = meanv),size=3)+
  geom_errorbar(aes(x = label,ymin = lower, ymax = upper),
                width = 0.6, # 控制上下两条短横线的长短
                size = 1 # 控制线条整体粗细
  )+theme(plot.title=element_text(hjust=0.5))+
  labs(x="m",y=TeX('mean value of \\mathbf{CV}(m)'),title =TeX('Setting 2,$\\Sigma_1$'))+
  xlim(0.5,10.5)+ylim(-0.3,2.2)
p4


err4eg3=read.csv("E:/first/eg3-44-cv.xls")[,-1]
err8eg3=read.csv("E:/first/eg3-88-cv-new.xls")[,-1]



df43 <- tibble(
  label = c(1:10),
  meanv = apply(err4eg3/100,2,mean),
  lower = meanv - apply(err4eg3/100,2,sd),
  upper = meanv + apply(err4eg3/100,2,sd),
  #group = c(rep("Group-1",7),rep("Group-2",7),rep("Group-3",8))
)


p5 <- ggplot(data = df43)+geom_point(aes(x = label, y = meanv),size=3)+
  geom_errorbar(aes(x = label,ymin = lower, ymax = upper),
                width = 0.6, # 控制上下两条短横线的长短
                size = 1 # 控制线条整体粗细
  )+theme(plot.title=element_text(hjust=0.5))+
  labs(x="m",y=TeX('mean value of \\mathbf{CV}(m)'),title =TeX('Setting 3,$\\Sigma_2$'))+
  xlim(0.5,10.5)+ylim(0,1.2)
p5


meanv=c()
lower=c()
upper=c()
for (i in c(1:10)) {
  meanv[i]=mean(err8eg3[which(abs(err8eg3[,i])<1000),i])/100
  lower[i]=meanv[i]-sd(err8eg3[which(abs(err8eg3[,i])<1000),i])/100
  upper[i]=meanv[i]+sd(err8eg3[which(abs(err8eg3[,i])<1000),i])/100
}


df83 <- tibble(
  label = c(1:10),
  meanv = meanv,
  lower = lower,
  upper = upper,
  #group = c(rep("Group-1",7),rep("Group-2",7),rep("Group-3",8))
)


p6 <- ggplot(data = df83)+geom_point(aes(x = label, y = meanv),size=3)+
  geom_errorbar(aes(x = label,ymin = lower, ymax = upper),
                width = 0.6, # 控制上下两条短横线的长短
                size = 1 # 控制线条整体粗细
  )+theme(plot.title=element_text(hjust=0.5))+
  labs(x="m",y=TeX('mean value of \\mathbf{CV}(m)'),title =TeX('Setting 3,$\\Sigma_1$'))+
  xlim(0.5,10.5)+ylim(-1.8,7.5)
p6




p7<-cowplot::plot_grid(p2,p1,p4,p3,p6,p5,nrow=3)
p7






