library(lme4)
library(lmerTest)

loadpath = '/Volumes/Tim/Photometry/10MfRatDataSet/'

triDf = read.csv(paste(loadpath,
                       'triframe4sameValtLeftChoiceReg.csv',sep=""))

altPathChooseLeftRegOut <- glmer( lrchoice ~ 1 + rt.1_left + rt.2_left + rt.3_left + rt.4_left + rt.5_left + 
                          (1 + rt.1_left + rt.2_left + rt.3_left + rt.4_left + rt.5_left | rat),
                          data = triDf[(triDf$samePath_t.1_left==0)&(triDf$leftLast==1),],family = binomial)

samePathChooseLeftRegOut <- glmer( lrchoice ~ 1 + rt.1_left + rt.2_left + rt.3_left + rt.4_left + rt.5_left + 
                                    (1 + rt.1_left + rt.2_left + rt.3_left + rt.4_left + rt.5_left | rat),
                                  data = triDf[(triDf$samePath_t.1_left==1)&(triDf$leftLast==1),],family = binomial)

samePathRegSum <- summary(samePathChooseLeftRegOut)
altPathRegSum <- summary(altPathChooseLeftRegOut)
samePathRegCoef = coef(samePathChooseLeftRegOut)
altPathRegCoef <- coef(altPathChooseLeftRegOut)

#write.csv(samePathDaRegCoef$`session:rat`,paste(
#  loadpath,"pChooseLeftSamePathCoefsBySessionAndRat.csv",sep=""))

write.csv(samePathRegCoef$rat,paste(loadpath,
                                      "pChooseLeftSamePathCoefsByRat.csv",sep=""))

#write.csv(altPathRegCoef$`session:rat`,paste(
#  loadpath,"pChooseLeftAltPathCoefsBySessionAndRat.csv",sep=""))

write.csv(altPathRegCoef$rat,paste(loadpath,
                                     "pChooseLeftAltPathCoefsByRat.csv",sep=""))

write.csv(samePathRegSum$coefficients,paste(
  loadpath,"pChooseLeftSamePathSum.csv",sep=""))

write.csv(altPathRegSum$coefficients,paste(
  loadpath,"pChooseLeftAltPathSum.csv",sep=""))