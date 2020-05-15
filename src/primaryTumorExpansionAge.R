setwd("dir/data")

library(ggplot2)
library(gmodels)
library(spatstat)


######colon
#time_step <- seq(-10,0,0.01)
#time_step <- seq(0,10,0.01)

pdf("Gompertz_growth_3cancers.pdf", width=4, height=4)
par(mar=c(5,5,3,2))

dat <- read.delim("DT_colon.txt", header=T)
dt_median <- weighted.median(dat$DT,dat$NumPatients/sum(dat$NumPatients))
print(dt_median)
dt_iqr <- weighted.quantile(dat$DT,dat$NumPatients/sum(dat$NumPatients), probs=seq(0,1,0.25), na.rm = TRUE)
dt_1qt <- as.numeric(dt_iqr[2])
dt_3qt <- as.numeric(dt_iqr[4])

limit_size <- 1e11
log_limit_size <- log(limit_size)
diameter <- 4.5
dx_size <- 4/3*pi*(diameter/2)^3*1e8

#for(dt in c(dt_1qt, dt_3qt, dt_median)){
for(dt in c(dt_median)) {
	beta <- 1/dt*log((log_limit_size-log(dx_size))/(log_limit_size-log(2*dx_size)))
	alpha <- log_limit_size*beta
	exp_age <- -1/beta*log(1-log(dx_size)*beta/alpha)
	exp_age <- exp_age/365
	print(c(beta,alpha))
	print(exp_age)
}


time_step <- seq(0,exp_age,0.01)
time_step2 <- time_step-exp_age
St <- exp(alpha/beta*(1-exp(-beta*time_step*365)))
plot(time_step2,log10(St),col="#3288bd",xlim=c(-8,0.5),ylim=c(0,12),type="l",lwd=3, yaxt="n", xlab="Years prior to diagnosis", ylab="Tumor size (# cells, log10)",cex.lab=1.3, cex.axis=1.2)
yxises <- c(0,2,4,6,8,10,12)
axis(2, at=yxises, labels=yxises, las=2, cex.axis=1.5)
abline(v=0, col="black",lwd=2,lty=2)
abline(v=time_step2[1], col="#3288bd",lty=2,lwd=3)
points(0,log10(dx_size),col="#3288bd",cex=1.5)

######breast
dat <- read.delim("DT_breast.txt", header=T)
dt_median <- weighted.median(dat$DT,dat$NumPatients/sum(dat$NumPatients))
print(dt_median)

dt_iqr <- weighted.quantile(dat$DT,dat$NumPatients/sum(dat$NumPatients), probs=seq(0,1,0.25), na.rm = TRUE)
dt_1qt <- as.numeric(dt_iqr[2])
dt_3qt <- as.numeric(dt_iqr[4])

limit_size <- 1e11
log_limit_size <- log(limit_size)
diameter <- 2.0
dx_size <- 4/3*pi*(diameter/2)^3*1e8

for(dt in c(dt_1qt, dt_3qt, dt_median)){
	beta <- 1/dt*log((log_limit_size-log(dx_size))/(log_limit_size-log(2*dx_size)))
	alpha <- log_limit_size*beta
	exp_age <- -1/beta*log(1-log(dx_size)*beta/alpha)
	exp_age <- exp_age/365
	print(c(beta,alpha))
	print(exp_age)
}


time_step <- seq(0,exp_age,0.01)
time_step2 <- time_step-exp_age

St <- exp(alpha/beta*(1-exp(-beta*time_step*365)))
lines(time_step2,log10(St),col="#9970ab",xlim=c(-8,0.5),ylim=c(0,12),lwd=3)
abline(v=time_step2[1], col="#9970ab",lty=2,lwd=3)
points(0,log10(dx_size),col="#9970ab",cex=1.5)

######lung
dat <- read.delim("DT_lung.txt", header=T)
dt_median <- weighted.median(dat$DT,dat$NumPatients/sum(dat$NumPatients))
print(dt_median)

dt_iqr <- weighted.quantile(dat$DT,dat$NumPatients/sum(dat$NumPatients), probs=seq(0,1,0.25), na.rm = TRUE)
dt_1qt <- as.numeric(dt_iqr[2])
dt_3qt <- as.numeric(dt_iqr[4])

limit_size <- 1e11
log_limit_size <- log(limit_size)
diameter <- 3.2
dx_size <- 4/3*pi*(diameter/2)^3*1e8

for(dt in c(dt_1qt, dt_3qt, dt_median)){
	beta <- 1/dt*log((log_limit_size-log(dx_size))/(log_limit_size-log(2*dx_size)))
	alpha <- log_limit_size*beta
	exp_age <- -1/beta*log(1-log(dx_size)*beta/alpha)
	exp_age <- exp_age/365
	print(c(beta,alpha))
	print(exp_age)
}

time_step <- seq(0,exp_age,0.01)
time_step2 <- time_step-exp_age

St <- exp(alpha/beta*(1-exp(-beta*time_step*365)))
lines(time_step2, log10(St),col="#d95f02",xlim=c(-8,0.5),ylim=c(0,12),lwd=3)
abline(v=time_step2[1], col="#d95f02",lty=2,lwd=3)
points(0,log10(dx_size),col="#d95f02",cex=1.5)

dev.off()

