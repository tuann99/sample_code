library(pwr)

# sample size calculation
mean_bh_time <- 54
sd_bh_time <- 23
effect_size <- 0.3 # moderate effect size
desired_power <- 0.8
alpha <- 0.05
sample_size <- pwr.r.test(
n = NULL, r = effect_size, sig.level = alpha, power = desired_power)
sample_size$n

# data
# 0 = male, 1 = female
# RRP = resting radial pulse in BPM (i.e., measurement via optical sensor for blood flow)
# RHR = resting heart rate in BPM (i.e., measurement via electrode)
# BH_time = breath hold time in seconds
# ID = subject number
path <- "C:\\Users\\nguye620\\OneDrive - Johns Hopkins\\School\\AS.410.645 (Biostatistics)\\final_project\\dat.csv"
dat <- read.csv(path, header = TRUE, row.names = 1)
head(dat)
summary(dat)

# number of men/women
men <- dat[dat$sex == 0,]
women <- dat[dat$sex == 1,]
nrow(men)
nrow(women)
mean(men$BH_time)
median(men$BH_time)
mean(women$BH_time)
median(women$BH_time)
mean(men$RRP)
median(men$RRP)
mean(women$RRP)
median(women$RRP)
mean(men$RHR)
median(men$RHR)
mean(women$RHR)
median(women$RHR)

# age groups
# dat$age <- as.integer(dat$age)
x <- dat[dat$age <= 22, ]
y <- dat[dat$age > 22, ]
nrow(x)
nrow(y)
mean(x$BH_time)
median(x$BH_time)
mean(y$BH_time)
median(y$BH_time)
mean(x$RRP)
median(x$RRP)
mean(y$RRP)
median(y$RRP)
mean(x$RHR)
median(x$RHR)
mean(y$RHR)
median(y$RHR)

# plot men vs women BH_times. Didn't really see a difference
# plot(men$RRP, men$BH_time, xlab = "RRP", ylab = "BH_time", main = "RRP vs BH_time for Men", col = "red")
# points(women$RRP, women$BH_time, col = "blue")
# legend("topright", legend=c('Men', 'Women'), col=c('red', 'blue'), pch = 15)

# plot RRPs vs BH_times
plot(dat$RRP, dat$BH_time, xlab = "Resting Radial Pulse (Beats Per Minute)", ylab = "Breath-hold Time (Seconds)", main = "Resting Radial Pulse vs Breath-Hold Time")

# plot RHRs vs BH_times
plot(dat$RHR, dat$BH_time, xlab = "Resting Heart Rate (Beats Per Minute)", ylab = "Breath-hold Time (Seconds)", main = "Resting Heart Rate vs Breath-Hold Time")

# plot RRPs and RHRs vs BH_times
# plot(dat$RRP, dat$BH_time,xlab = "Heart Rate (BPM)", ylab = "BH Time (Seconds)", main = "RRP and RHR vs. BH Time", col = "red")
# points(dat$RHR, dat$BH_time, col = "blue")
# legend("topright", legend=c('RRP', 'RHR'), col=c('red', 'blue'), pch = 15)

# check for outliers via IQR method
# boxplot(dat$RRP, main = "RRP", ylab = "RRP")
# boxplot(dat$RHR, main = "RHR", ylab = "RHR")
# RRP outliers
quartile_1_RRP<- quantile(dat$RRP, 0.25)
quartile_3_RRP <- quantile(dat$RRP, 0.75)
iqr_RRP <- quartile_3_RRP - quartile_1_RRP
lower_bound_RRP <- quartile_1_RRP - 1.5 * iqr_RRP
upper_bound_RRP <- quartile_3_RRP + 1.5 * iqr_RRP
outliers_RRP <- dat$RRP[dat$RRP < lower_bound_RRP | dat$RRP > upper_bound_RRP]

# RHR outliers
quartile_1_RHR <- quantile(dat$RHR, 0.25)
quartile_3_RHR <- quantile(dat$RHR, 0.75)
iqr_RHR <- quartile_3_RHR - quartile_1_RHR
lower_bound_RHR <- quartile_1_RHR - 1.5 * iqr_RHR
upper_bound_RHR <- quartile_3_RHR + 1.5 * iqr_RHR
outliers_RHR <- dat$RHR[dat$RHR < lower_bound_RHR | dat$RHR > upper_bound_RHR]

outliers_RRP
outliers_RHR

# # analysis via pearson correlation hypothesis test
# H0: there is no correlation between BH_time and RHR
# HA: there is a negative correlation between BH_time and RHR
rhr_cor <- cor.test(dat$BH_time, dat$RHR, method = "pearson", alternative = "less")
rrp_cor <- cor.test(dat$BH_time, dat$RRP, method = "pearson", alternative = "less")
rhr_cor
rrp_cor

# analysis via linear regression
rhr_reg <- lm(dat$BH_time ~ dat$RHR)
rrp_reg <- lm(dat$BH_time ~ dat$RRP)
summary(rhr_reg)
summary(rrp_reg)

# lm adjusting for age and sex
rrp_reg_cov <- lm(dat$BH_time ~ dat$RRP + dat$sex + dat$age)
summary(rrp_reg_cov)

rhr_reg_cov <- lm(dat$BH_time ~ dat$RHR + dat$sex + dat$age)
summary(rhr_reg_cov)
