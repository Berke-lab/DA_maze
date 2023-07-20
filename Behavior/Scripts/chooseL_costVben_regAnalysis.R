library(lme4)
library(lmerTest)
library(scales)

loadpath = "/Volumes/Tim/Photometry/10MfRatDataSet/"

choiceDf = read.csv(paste(loadpath,'costVbenDf4reg.csv',sep=""))

choiceDf$pRwdScaled <- rescale(choiceDf$pRwdDif,c(0,1))
choiceDf$ldifScaled <- rescale(choiceDf$ldif,c(0,1))

chooseLeftRegOut <- glmer(choose_L ~ 1 + pRwdScaled + ldifScaled + 
                    (1 + pRwdScaled + ldifScaled |rat)+
                      (1 + pRwdScaled + ldifScaled | session:rat),
                  data = choiceDf[choiceDf$tri>25,],family = binomial)

chooseLcoefs <- coef(chooseLeftRegOut)
chooseLsum <- summary(chooseLeftRegOut)
write.csv(chooseLsum$coefficients,paste(loadpath,"chooseLregSummary.csv",sep=''))
write.csv(chooseLcoefs$rat,paste(loadpath,"chooseLregCoefsByRat.csv",sep=''))