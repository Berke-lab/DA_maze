library(lme4)
library(lmerTest)
library(scales)

loadpath = "/Volumes/Tim/Photometry/10MfRatDataSet/"

daOptGamValDf <- read.csv(paste(loadpath,"valDaDf.csv",sep=""))

daOptGamValDf$vel_scaled <- rescale(daOptGamValDf$vel,c(0,1))

valRegOut <- lmer(DA ~ 1 + Value + vel_scaled + (1 + Value + vel_scaled | rat) + 
                 (1 + Value + vel_scaled | session:rat),data = daOptGamValDf)

DaRegSum <- summary(valRegOut)
DaRegCoef = coef(valRegOut)
write.csv(DaRegCoef$`session:rat`,paste(loadpath,"valDACoefsBySessionAndRat.csv",sep=""))
write.csv(DaRegCoef$rat,paste(loadpath,"valDACoefsByRat.csv",sep=""))
write.csv(DaRegSum$coefficients,paste(loadpath,"valDASum.csv",sep=""))

#compute nested regressions

velRegOut <- lmer(DA ~ 1 + vel_scaled + (1 + vel_scaled | rat) + 
                    (1 + vel_scaled | session:rat),data = daOptGamValDf)

valOnlyRegOut <- lmer(DA ~ 1 + Value + (1 + Value | rat) + 
                        (1 + Value | session:rat),data = daOptGamValDf)

#compute coefficients of partial determination

sse_valAndVel <- sum((fitted(valRegOut) - daOptGamValDf$DA)^2)

sse_vel <- sum((fitted(velRegOut) - daOptGamValDf$DA)^2)

sse_val <- sum((fitted(valOnlyRegOut) - daOptGamValDf$DA)^2)

ssr_valAndVel <- sum((fitted(valRegOut) - mean(daOptGamValDf$DA))^2)

ssr_val <- sum((fitted(valOnlyRegOut) - mean(daOptGamValDf$DA))^2)

ssr_vel <- sum((fitted(velRegOut) - mean(daOptGamValDf$DA))^2)

valPartialR2 <- (ssr_valAndVel - ssr_vel)/sse_vel

velPartialR2 <- (ssr_valAndVel - ssr_val)/sse_val

partialR2df <- data.frame (value  = c(valPartialR2),
                  velocity = c(velPartialR2))

write.csv(partialR2df,paste(loadpath,"valDaPartialR2s.csv",sep=""))