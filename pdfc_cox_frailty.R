
dt <- read.csv("dwell_time.csv")
set.seed(0)
list.of.packages <- c("survival")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library("survival")
  
intervals <- read.csv("DTintervals.csv")

## frailty model
# compute intervals from start time and end time
dt$int = dt$end_time-dt$start_time
eFrail.coxph <- coxph(Surv(int,event) ~ cov1 + frailty(Subject), data = dt)

columnnames = c("time")
for (c in 1:length(intervals$Var1)){
  columnnames <- append(columnnames,paste("Val",toString(c),sep="_"))
}

## ** expected number of transition per subject
lambda0 <- basehaz(eFrail.coxph, centered = TRUE)$hazard
eXbeta <- predict(eFrail.coxph, type = "risk", newdata = data.frame(cov1=intervals$Var1))
z <- eFrail.coxph$frail

cumhazardFrail <- tcrossprod(lambda0,eXbeta) %o% exp(z)
Mhaz.exp <- cbind(basehaz(eFrail.coxph, centered = TRUE)$time,apply(cumhazardFrail,1:2,mean))
colnames(Mhaz.exp) <- columnnames

survFrail <- exp(-cumhazardFrail)
Msurv.exp <- cbind(basehaz(eFrail.coxph, centered = TRUE)$time,apply(survFrail,1:2,mean))
colnames(Msurv.exp) <- columnnames

Msurv.exp2 = as.data.frame(Msurv.exp)
write.csv(Msurv.exp2,file="dwell_time_Surv_curves.csv",na="NaN")

survsum = summary(eFrail.coxph)
beta <- survsum$conf.int[1]
betaCI1 <- survsum$conf.int[3]
betaCI2 <- survsum$conf.int[4]
nvar <- length(beta)
nfrail <- nrow(eFrail.coxph$var) - nvar
se <- survsum$coefficients[1,2]
z<- round(beta/se, 2)
p <- survsum$coefficients[1,6]
a_df=data.frame(cbind(beta,betaCI1,betaCI2,se,z,p))

write.csv(a_df,file="dwell_time_Surv_stats.csv",na="NaN")


