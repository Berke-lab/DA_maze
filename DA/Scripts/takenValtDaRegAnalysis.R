library(lme4)
library(lmerTest)

loadpath = "/Volumes/Tim/Photometry/10MfRatDataSet/"

hexDfTakenVAlt <- read.csv(paste(loadpath,'hexDf_with_PreCp_and_prevRwdAltSame_correctAlignment.csv',sep=""))
preCpSame2Df <- hexDfTakenVAlt[(hexDfTakenVAlt$preCP==1)&(hexDfTakenVAlt$same2==1),]
preCpDaAltPathRegOut <- lmer( DA ~ 1 + rt.1 + rt.2 + rt.3 + rt.4 + rt.5 + 
                          (1 + rt.1 + rt.2 + rt.3 + rt.4 + rt.5 | rat),
                        data = preCpSame2Df[preCpSame2Df$samePath_t.1==0,])

preCpDaSamePathRegOut <- lmer( DA ~ 1 + rt.1 + rt.2 + rt.3 + rt.4 + rt.5 + 
                          (1 + rt.1 + rt.2 + rt.3 + rt.4 + rt.5 | rat),
                          data = preCpSame2Df[preCpSame2Df$samePath_t.1==1,])

samePathDaRegSum <- summary(preCpDaSamePathRegOut)
altPathDaRegSum <- summary(preCpDaAltPathRegOut)
samePathDaRegCoef = coef(preCpDaSamePathRegOut)
altPathDaRegCoef <- coef(preCpDaAltPathRegOut)
write.csv(samePathDaRegCoef$`session:rat`,paste(loadpath,"DApreCpSamePathCoefsBySessionAndRat.csv",sep=""))
write.csv(samePathDaRegCoef$rat,paste(loadpath,"DApreCpSamePathCoefsByRat.csv",sep=""))
write.csv(altPathDaRegCoef$`session:rat`,paste(loadpath,"DApreCpAltPathCoefsBySessionAndRat.csv",sep=""))
write.csv(altPathDaRegCoef$rat,paste(loadpath,"DApreCpAltPathCoefsByRat.csv",sep=""))
write.csv(samePathDaRegSum$coefficients,paste(loadpath,"DApreCpSamePathSum.csv",sep=""))
write.csv(altPathDaRegSum$coefficients,paste(loadpath,"DApreCpAltPathSum.csv",sep=""))