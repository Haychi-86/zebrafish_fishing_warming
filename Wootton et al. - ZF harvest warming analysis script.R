# This code is associated with the publication "Harvesting that preserves the large fish might mitigate the worst impacts of warming".
# It analyses life history data recorded during a multigenerational experiment on zebrafish that investigated the impacts of fisheries harvests and climate warming. 

# see the manuscript for a full description and interpretation of the analyses. Associated data with this analysis contains the variables 'Deviance' (calculated as 
# the difference between each record and the mean of the relevant control in each generation), 'Temperature' ('L' for 26 °C and 'H' for 30 °C treatments), 'Selection' ('C' 
# for uniform size selection controls, 'G' for Gaussian fisheries selection and 'S' for sigmoid fisheries selection), Tank (numbered 1-6 indicating housing tanks), 'Population'
# (specific codes indicating differing replicate populations) and 'Generation' (numbered 0-7 for generations PG2 to CG2). A factor for specific combinations of fisheries selection and warming treatments is 
# generated in the code below for each analysis.


#### Start analysis ####

rm(list = ls()) # clear R environment
#setwd("C:Users/...") # Set the directory to your location of choice
# Or, set the working directory to the location of this file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

dir.create("InputsZf")  # Create directory structure - place data files in this 'inputs' folder
library(lme4)
library(MQMF)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(effects)

#### Adult length analysis ####

# Compare adult total lengths across treatments and generations. 

# load and prepare data for analysis during 'F' generations
length_data_dev <- read.csv(file = "InputsZf/Length_Dev_data.csv")

Dev_length_dataF <- length_data_dev %>% filter (Generation > -1 & Generation < 5 )

Dev_length_dataF$Treat <- as.factor(paste(Dev_length_dataF$Temperature,Dev_length_dataF$Selection, sep = "_"))

Dev_length_dataF$Treat <- relevel(Dev_length_dataF$Treat , ref="L_C")

# first fit the full model to check fit
LengthFM <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataF, REML = T)# 

summary(LengthFM)

plot(LengthFM) # check resid spread 

qqnorm(residuals(LengthFM))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=Dev_length_dataF$Generation,pearson=residuals(LengthFM,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Treat=Dev_length_dataF$Treat,pearson=residuals(LengthFM,type="pearson")),aes(x=Treat,y=pearson)) + geom_point() + theme_bw()

# now run the model selection for Deviance. 
# specify the full model
LengthFMa <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataF, REML = F) 
# specify the additive model
LengthFMb <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataF, REML = F) 

# test the models using anova
anova(LengthFMa, LengthFMb)
#Data: Dev_length_dataF
#Models:
#LengthFMb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#LengthFMa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#          npar   AIC   BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#LengthFMb   11 14748 14812 -7362.8    14726                         
#LengthFMa   16 14717 14812 -7342.6    14685 40.434  5  1.221e-07 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
LengthFMbest <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataF, REML = T) 

plot(allEffects(LengthFMbest)) # gives good visual of parameter estimates

summary(LengthFMbest)
#Fixed effects:
#                      Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)          1.402e-11  9.114e-01  7.216e+00   0.000 1.000000    
#TreatH_C            -1.880e+00  8.196e-01  1.245e+01  -2.294 0.039926 *  
#TreatH_G            -2.022e+00  8.196e-01  1.245e+01  -2.467 0.028999 *  
#TreatH_S            -3.827e+00  8.196e-01  1.245e+01  -4.669 0.000492 ***
#TreatL_G            -1.111e-01  6.467e-01  1.553e+01  -0.172 0.865808    
#TreatL_S            -2.304e+00  6.467e-01  1.553e+01  -3.563 0.002700 ** 
#Generation          -3.435e-12  3.118e-01  3.932e+00   0.000 1.000000    
#TreatH_C:Generation  2.511e-01  1.718e-01  2.673e+03   1.461 0.144031    
#TreatH_G:Generation  7.856e-01  1.718e-01  2.673e+03   4.572 5.06e-06 ***
#TreatH_S:Generation  1.589e-01  1.718e-01  2.673e+03   0.925 0.355220    
#TreatL_G:Generation  4.944e-01  1.718e-01  2.673e+03   2.877 0.004040 ** 
#TreatL_S:Generation -1.778e-01  1.718e-01  2.673e+03  -1.035 0.300950        

confint(LengthFMbest)
#                          2.5 %     97.5 %
#.sig01               0.24409305  0.8139911
#.sig02               0.05299414  1.2439779
#.sig03               0.41013866  1.6180261
#.sigma               3.54634055  3.7414874
#(Intercept)         -1.62959544  1.6296023
#TreatH_C            -3.37404815 -0.3859457
#TreatH_G            -3.51627038 -0.5281679
#TreatH_S            -5.32071482 -2.3326124
#TreatL_G            -1.25633378  1.0341120
#TreatL_S            -3.44966712 -1.1592213
#Generation          -0.59469202  0.5946884
#TreatH_C:Generation -0.08548270  0.5877049
#TreatH_G:Generation  0.44896174  1.1221494
#TreatH_S:Generation -0.17770492  0.4954827
#TreatL_G:Generation  0.15785063  0.8310383
#TreatL_S:Generation -0.51437159  0.1588160

# Lets look at F5 generation to see the effects after treatments were relaxed

# filter and prepare data for F5 generation analysis
Dev_length_dataF5 <- length_data_dev %>% filter (Generation == 5)

Dev_length_dataF5$Treat <- as.factor(paste(Dev_length_dataF5$Temperature,Dev_length_dataF5$Selection, sep = "_"))

Dev_length_dataF5$Treat <- relevel(Dev_length_dataF5$Treat , ref="L_C")

# first fit the full model to check fit
LengthF5M <- lmer(Deviance ~ Treat   + (1|Tank), data = Dev_length_dataF5, REML = T)

summary(LengthF5M)

plot(LengthF5M) # check resid spread 

qqnorm(residuals(LengthF5M))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Treat=Dev_length_dataF5$Treat,pearson=residuals(LengthF5M,type="pearson")),aes(x=Treat,y=pearson)) + geom_point() + theme_bw()

# now run the model selection for Deviance. 
# specify the full model
LengthF5Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = Dev_length_dataF5, REML = F) 
# specify the null model
LengthF5Mb <- lmer(Deviance ~    + (1|Tank), data = Dev_length_dataF5, REML = F) 

# test the models using anova
anova(LengthF5Ma, LengthF5Mb)
#Data: Dev_length_dataF5
#Models:
#LengthF5Mb: Deviance ~ +(1 | Tank)
#LengthF5Ma: Deviance ~ Treat + (1 | Tank)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#LengthF5Mb    3 2735.2 2747.8 -1364.6   2729.2                         
#LengthF5Ma    8 2609.6 2643.4 -1296.8   2593.6 135.58  5  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
LengthF5Mbest <- lmer(Deviance ~ Treat  + (1|Tank) , data = Dev_length_dataF5, REML = T) 

plot(allEffects(LengthF5Mbest)) # gives good visual of parameter estimates

summary(LengthF5Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)  2.535e-15  4.521e-01  8.746e+00   0.000  1.00000    
#TreatH_C    -1.344e+00  6.394e-01  8.746e+00  -2.103  0.06571 .  
#TreatH_G    -2.000e-01  6.394e-01  8.746e+00  -0.313  0.76175    
#TreatH_S    -3.021e+00  6.917e-01  1.127e+01  -4.368  0.00106 ** 
#TreatL_G     3.722e+00  4.580e-01  5.000e+02   8.127 3.48e-15 ***
#TreatL_S    -8.222e-01  4.580e-01  5.000e+02  -1.795  0.07323 .      

confint(LengthF5Mbest)
#.                2.5 %      97.5 %
#.sig01       0.0000000  0.99083947
#.sigma       2.8807954  3.25957157
#(Intercept) -0.8312368  0.83123710
#TreatH_C    -2.5199908 -0.16889766
#TreatH_G    -1.3755463  0.97554678
#TreatH_S    -4.2948116 -1.74534074
#TreatL_G     2.8263327  4.61811175
#TreatL_S    -1.7181118  0.07366731

# CG1 analysis

# filter and prepare data for CG1 generation analysis
Dev_length_dataCG1 <- length_data_dev %>% filter (Generation == 6 )

Dev_length_dataCG1$Treat <- as.factor(paste(Dev_length_dataCG1$Temperature,Dev_length_dataCG1$Selection, sep = "_"))

Dev_length_dataCG1$Treat <- relevel(Dev_length_dataCG1$Treat , ref="L_C")

# first fit the full model to check fit
LengthCG1M <- lmer(Deviance ~ Treat  + (1|Tank), data = Dev_length_dataCG1, REML = T)

summary(LengthCG1M)

plot(LengthCG1M) # check resid spread 

qqnorm(residuals(LengthCG1M))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Treat=Dev_length_dataCG1$Treat,pearson=residuals(LengthCG1M,type="pearson")),aes(x=Treat,y=pearson)) + geom_point() + theme_bw()

# now run the model selection for Deviance. 
# specify the full model
LengthCG1Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = Dev_length_dataCG1, REML = F) 
# specify the null model
LengthCG1Mb <- lmer(Deviance ~    + (1|Tank), data = Dev_length_dataCG1, REML = F) 

# test the models using anova
anova(LengthCG1Ma, LengthCG1Mb)
#Data: Dev_length_dataCG1
#Models:
#LengthCG1Mb: Deviance ~ +(1 | Tank)
#LengthCG1Ma: Deviance ~ Treat + (1 | Tank)
#            npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#LengthCG1Mb    3 2630.5 2643.2 -1312.2   2624.5                         
#LengthCG1Ma    8 2551.8 2585.6 -1267.9   2535.8 88.714  5  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
LengthCG1Mbest <- lmer(Deviance ~ Treat  + (1|Tank) , data = Dev_length_dataCG1, REML = T) 

plot(allEffects(LengthCG1Mbest)) # gives good visual of parameter estimates

summary(LengthCG1Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept) -4.311e-14  3.859e-01  1.170e+01   0.000 1.000000    
#TreatH_C    -2.244e+00  5.458e-01  1.170e+01  -4.112 0.001517 ** 
#TreatH_G    -2.444e-01  5.458e-01  1.170e+01  -0.448 0.662430    
#TreatH_S    -2.582e+00  5.980e-01  1.515e+01  -4.318 0.000596 ***
#TreatL_G     2.467e+00  4.336e-01  5.002e+02   5.689 2.18e-08 ***
#TreatL_S    -4.000e-01  4.336e-01  5.002e+02  -0.922 0.356715       

confint(LengthCG1Mbest)
#.                2.5 %     97.5 %
#.sig01       0.0000000  0.7754434
#.sigma       2.7270484  3.0854540
#(Intercept) -0.7107699  0.7107702
#TreatH_C    -3.2496248 -1.2392636
#TreatH_G    -1.2496248  0.7607364
#TreatH_S    -3.6836514 -1.4771483
#TreatL_G     1.6186118  3.3147216
#TreatL_S    -1.2480549  0.4480549

# CG2 analysis

# filter and prepare data for CG2 generation analysis
Dev_length_dataCG2 <- length_data_dev %>% filter (Generation > 6 )

Dev_length_dataCG2$Treat <- as.factor(paste(Dev_length_dataCG2$Temperature,Dev_length_dataCG2$Selection, sep = "_"))

Dev_length_dataCG2$Treat <- relevel(Dev_length_dataCG2$Treat , ref="L_C")

# first fit the full model to check fit
LengthCG2M <- lmer(Deviance ~ Treat  + (1|Tank), data = Dev_length_dataCG2, REML = T)

summary(LengthCG2M)

plot(LengthCG2M) # check resid spread 

qqnorm(residuals(LengthCG2M))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Treat=Dev_length_dataCG2$Treat,pearson=residuals(LengthCG2M,type="pearson")),aes(x=Treat,y=pearson)) + geom_point() + theme_bw()

# now run the model selection for Deviance. 
# specify the full model
LengthCG2Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = Dev_length_dataCG2, REML = F) 
# specify the null model
LengthCG2Mb <- lmer(Deviance ~    + (1|Tank), data = Dev_length_dataCG2, REML = F) 

# test the models using anova
anova(LengthCG2Ma, LengthCG2Mb)
#Data: Dev_length_dataCG2
#Models:
#LengthCG2Mb: Deviance ~ +(1 | Tank)
#LengthCG2Ma: Deviance ~ Treat + (1 | Tank)
#            npar    AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)    
#LengthCG2Mb    3 2588.0 2600.7  -1291   2582.0                         
#LengthCG2Ma    8 2507.9 2541.8  -1246   2491.9 90.141  5  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
LengthCG2Mbest <- lmer(Deviance ~ Treat   + (1|Tank) , data = Dev_length_dataCG2, REML = T) 

plot(allEffects(LengthCG2Mbest)) # gives good visual of parameter estimates

summary(LengthCG2Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept) -3.579e-14  4.267e-01  8.448e+00   0.000  1.00000    
#TreatH_C    -8.444e-01  6.034e-01  8.448e+00  -1.399  0.19730    
#TreatH_G     5.444e-01  6.034e-01  8.448e+00   0.902  0.39190    
#TreatH_S    -1.590e+00  6.495e-01  1.080e+01  -2.449  0.03269 *  
#TreatL_G     2.233e+00  4.142e-01  5.002e+02   5.391 1.08e-07 ***
#TreatL_S    -1.233e+00  4.142e-01  5.002e+02  -2.977  0.00305 **        

confint(LengthCG2Mbest)
#                  2.5 %     97.5 %
#.sig01       0.09968492  0.9587051
#.sigma       2.60530398  2.9477341
#(Intercept) -0.78799975  0.7880001
#TreatH_C    -1.95884438  0.2699559
#TreatH_G    -0.56995549  1.6588448
#TreatH_S    -2.78792360 -0.3934290
#TreatL_G     1.42313497  3.0435317
#TreatL_S    -2.04353170 -0.4231350

# now compare adult total lengths across treatments and generations in warmed populations only.
# this analysis compares fisheries selection treatments to fisheries size selection controls (uniform size selection) in warmed populations only

# load and prepare data for warmed 'F' generation analysis
Dev_length_dataFH <- read.csv(file = "InputsZf/Length_Dev_data_Warmed.csv")

Dev_length_dataFH$Treat <- as.factor(paste(Dev_length_dataFH$Temperature,Dev_length_dataFH$Selection, sep = "_"))

Dev_length_dataFH$Treat <- relevel(Dev_length_dataFH$Treat , ref="H_C")

# first fit the full model to check fit
LengthFMH <- lmer(Deviance ~ Treat*Generation   + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataFH, REML = T)# 

summary(LengthFMH)

plot(LengthFMH) # check resid spread 

qqnorm(residuals(LengthFMH))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=Dev_length_dataFH$Generation,pearson=residuals(LengthFMH,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Treat=Dev_length_dataFH$Treat,pearson=residuals(LengthFMH,type="pearson")),aes(x=Treat,y=pearson)) + geom_point() + theme_bw()

# now run the model selection for Deviance. 
# specify the full model
LengthFMHa <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataFH, REML = F) 
# specify the additive model
LengthFMHb <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataFH, REML = F) 

# test the models using anova
anova(LengthFMHa, LengthFMHb)
#Data: Dev_length_dataFH
#Models:
#LengthFMHb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#LengthFMHa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#LengthFMHb    8 7162.7 7204.4 -3573.4   7146.7                         
#LengthFMHa   10 7148.6 7200.7 -3564.3   7128.6 18.073  2   0.000119 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
LengthFMHbest <- lmer(Deviance ~ Treat*Generation   + (1|Population) + (1|Tank) + (1|Generation), data = Dev_length_dataFH, REML = T)# 

plot(allEffects(LengthFMHbest)) # gives good visual of parameter estimates

summary(LengthFMHbest)
#Fixed effects:
#                      Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)         -3.772e-13  7.056e-01  6.084e+00   0.000 1.000000    
#TreatH_G            -1.422e-01  5.523e-01  8.917e+00  -0.258 0.802634    
#TreatH_S            -1.947e+00  5.523e-01  8.917e+00  -3.525 0.006564 ** 
#Generation           5.793e-14  1.837e-01  5.316e+00   0.000 1.000000    
#TreatH_G:Generation  5.344e-01  1.588e-01  1.335e+03   3.367 0.000783 ***
#TreatH_S:Generation -9.222e-02  1.588e-01  1.335e+03  -0.581 0.561389      

confint(LengthFMHbest)
#                           2.5 %     97.5 %
#.sig01               0.05469044  0.8618527
#.sig02               0.13631605  0.9626354
#.sig03               0.17739841  2.1365211
#.sigma               3.24150804  3.4968340
#(Intercept)         -1.35574435  1.3557411
#TreatH_G            -1.16301496  0.8785695
#TreatH_S            -2.96745940 -0.9258750
#Generation          -0.36684296  0.3668435
#TreatH_G:Generation  0.22330370  0.8455852
#TreatH_S:Generation -0.40336296  0.2189185

#### Adult weight analysis ####

# Compare adult weights across treatments and generations. 

# load and prepare data for analysis during 'F' generations
Weight_data_dev <- read.csv(file = "InputsZf/Weight_Dev_data.csv")

Weight_data_devF <- Weight_data_dev %>% filter (Generation > -1 & Generation < 5 )

Weight_data_devF$Treat <- as.factor(paste(Weight_data_devF$Temperature,Weight_data_devF$Selection, sep = "_"))

Weight_data_devF$Treat <- relevel(Weight_data_devF$Treat , ref="L_C")

# first fit the full model to check fit
WeightFM <- lmer(Deviance ~ Treat*Generation   + (1|Population) + (1|Tank) + (1|Generation), data = Weight_data_devF, REML = T)# 

summary(WeightFM)

plot(WeightFM) # check resid spread 

qqnorm(residuals(WeightFM))

# now run the model selection for Deviance. 
# specify the full model
WeightFMa <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Weight_data_devF, REML = F) 
# specify the additive model
WeightFMb <- lmer(Deviance ~ Treat + Generation + (1|Population) + (1|Tank) + (1|Generation), data = Weight_data_devF, REML = F) 

# test the models using anova
anova(WeightFMa, WeightFMb)
#Data: Weight_data_devF
#Models:
#WeightFMb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#WeightFMa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#WeightFMb   11 1725.8 1790.7 -851.90   1703.8                         
#WeightFMa   16 1703.1 1797.5 -835.57   1671.1 32.665  5  4.385e-06 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
WeightFMbest <- lmer(Deviance ~ Treat*Generation   + (1|Population) + (1|Tank) + (1|Generation), data = Weight_data_devF, REML = T)# 

plot(allEffects(WeightFMbest)) # gives good visual of parameter estimates

summary(WeightFMbest)
#Fixed effects:
#                      Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)          3.549e-14  8.605e-02  5.588e+00   0.000  1.00000    
#TreatH_C            -1.980e-01  6.594e-02  1.732e+01  -3.002  0.00789 ** 
#TreatH_G            -1.657e-01  6.594e-02  1.732e+01  -2.513  0.02212 *  
#TreatH_S            -3.094e-01  6.594e-02  1.732e+01  -4.693  0.00020 ***
#TreatL_G             1.018e-02  5.885e-02  1.522e+01   0.173  0.86494    
#TreatL_S            -1.375e-01  5.885e-02  1.522e+01  -2.336  0.03357 *  
#Generation          -2.412e-14  3.148e-02  3.705e+00   0.000  1.00000    
#TreatH_C:Generation  7.982e-03  1.544e-02  2.673e+03   0.517  0.60512    
#TreatH_G:Generation  3.907e-02  1.544e-02  2.673e+03   2.531  0.01143 *  
#TreatH_S:Generation -1.159e-02  1.544e-02  2.673e+03  -0.751  0.45264    
#TreatL_G:Generation  4.790e-02  1.544e-02  2.673e+03   3.103  0.00193 ** 
#TreatL_S:Generation -2.250e-02  1.544e-02  2.673e+03  -1.458  0.14502   

confint(WeightFMbest)
#                           2.5 %       97.5 %
#.sig01               0.022897124  0.074699547
#.sig02               0.000000000  0.087587784
#.sig03               0.042221499  0.162504335
#.sigma               0.318573641  0.336104002
#(Intercept)         -0.154769120  0.154768957
#TreatH_C            -0.316820881 -0.079138867
#TreatH_G            -0.284556420 -0.046874422
#TreatH_S            -0.428271976 -0.190589978
#TreatL_G            -0.093962736  0.114322728
#TreatL_S            -0.241604959 -0.033319494
#Generation          -0.059203539  0.059203161
#TreatH_C:Generation -0.022254556  0.038219001
#TreatH_G:Generation  0.008832111  0.069305667
#TreatH_S:Generation -0.041831223  0.018642334
#TreatL_G:Generation  0.017664333  0.078137889
#TreatL_S:Generation -0.052739001  0.007734556

# Lets look at F5 generation to see the effects after treatments were relaxed

# filter and prepare data for F5 generation analysis
Weight_data_devF5 <- Weight_data_dev %>% filter (Generation == 5 )

Weight_data_devF5$Treat <- as.factor(paste(Weight_data_devF5$Temperature,Weight_data_devF5$Selection, sep = "_"))

Weight_data_devF5$Treat <- relevel(Weight_data_devF5$Treat , ref="L_C")

# first fit the full model to check fit
WeightF5M <- lmer(Deviance ~ Treat + (1|Tank), data = Weight_data_devF5, REML = T)

summary(WeightF5M)

plot(WeightF5M) # check resid spread 

qqnorm(residuals(WeightF5M))

# now run the model selection for Deviance. 
# specify the full model
WeightF5Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = Weight_data_devF5, REML = F) 
# specify the null model
WeightF5Mb <- lmer(Deviance ~    + (1|Tank), data = Weight_data_devF5, REML = F) 

# test the models using anova
anova(WeightF5Ma, WeightF5Mb)
#Data: Weight_data_devF5
#Models:
#WeightF5Mb: Deviance ~ +(1 | Tank)
#WeightF5Ma: Deviance ~ Treat + (1 | Tank)
#           npar    AIC    BIC   logLik deviance  Chisq Df Pr(>Chisq)    
#WeightF5Mb    3 302.72 315.42 -148.359   296.72                         
#WeightF5Ma    8 190.21 224.08  -87.103   174.21 122.51  5  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
Wgt_F5Mbest <- lmer(Deviance ~  Treat  + (1|Tank) , data = Weight_data_devF5, REML = T) 

plot(allEffects(Wgt_F5Mbest)) # gives good visual of parameter estimates

summary(Wgt_F5Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)  2.415e-16  3.894e-02  1.007e+01   0.000 1.000000    
#TreatH_C    -1.710e-01  5.506e-02  1.007e+01  -3.105 0.011061 *  
#TreatH_G    -1.139e-01  5.506e-02  1.007e+01  -2.069 0.065210 .  
#TreatH_S    -2.807e-01  6.016e-02  1.306e+01  -4.666 0.000436 ***
#TreatL_G     3.064e-01  4.280e-02  4.999e+02   7.160 2.91e-12 ***
#TreatL_S    -1.287e-01  4.280e-02  4.999e+02  -3.008 0.002765 ** 

confint(Wgt_F5Mbest)
#                  2.5 %      97.5 %
#.sig01       0.00000000  0.08018993
#.sigma       0.26921583  0.30462589
#(Intercept) -0.07122932  0.07122935
#TreatH_C    -0.27172236 -0.07025538
#TreatH_G    -0.21465569 -0.01318871
#TreatH_S    -0.39108875 -0.16876258
#TreatL_G     0.22269789  0.39014655
#TreatL_S    -0.21244655 -0.04499789

# CG2 analysis
# note weights were not recorded in CG1 due to restrictions around coronavirus lockdowns. See the manuscript for more detail

# filter and prepare data for CG2 generation analysis
Weight_data_devCG2 <- Weight_data_dev %>% filter (Generation > 6 )

Weight_data_devCG2$Treat <- as.factor(paste(Weight_data_devCG2$Temperature,Weight_data_devCG2$Selection, sep = "_"))

Weight_data_devCG2$Treat <- relevel(Weight_data_devCG2$Treat , ref="L_C")

# fit the full model
WeightCG2M <- lmer(Deviance ~ Treat + (1|Tank), data = Weight_data_devCG2, REML = T)

summary(WeightCG2M)

plot(WeightCG2M) # check resid spread 

qqnorm(residuals(WeightCG2M))

# now run the model selection for Deviance. 
# specify the full model
WeightCG2Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = Weight_data_devCG2, REML = F) 
# specify the null model
WeightCG2Mb <- lmer(Deviance ~     + (1|Tank), data = Weight_data_devCG2, REML = F) 

# test the models using anova
anova(WeightCG2Ma,WeightCG2Mb)
#Data: Weight_data_devCG2
#Models:
#WeightCG2Mb: Deviance ~ +(1 | Tank)
#WeightCG2Ma: Deviance ~ Treat + (1 | Tank)
#            npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#WeightCG2Mb    3 293.25 305.95 -143.62   287.25                         
#WeightCG2Ma    8 233.05 266.92 -108.52   217.05 70.204  5  9.294e-14 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
WeightCG2Mbest <- lmer(Deviance ~  Treat  + (1|Tank) , data = Weight_data_devCG2, REML = T) 

plot(allEffects(WeightCG2Mbest)) # gives good visual of parameter estimates

summary(WeightCG2Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)  4.803e-16  3.858e-02  1.257e+01   0.000 1.000000    
#TreatH_C    -1.259e-01  5.457e-02  1.257e+01  -2.307 0.038771 *  
#TreatH_G    -3.294e-02  5.457e-02  1.257e+01  -0.604 0.556738    
#TreatH_S    -1.750e-01  6.000e-02  1.620e+01  -2.917 0.009975 ** 
#TreatL_G     1.841e-01  4.468e-02  5.002e+02   4.120 4.43e-05 ***
#TreatL_S    -1.585e-01  4.468e-02  5.002e+02  -3.548 0.000424 ***

confint(WeightCG2Mbest)
#                  2.5 %      97.5 %
#.sig01       0.00000000  0.07516634
#.sigma       0.28100408  0.31793917
#(Intercept) -0.07097450  0.07097453
#TreatH_C    -0.22626199 -0.02551574
#TreatH_G    -0.13331755  0.06742870
#TreatH_S    -0.28540344 -0.06371508
#TreatL_G     0.09669087  0.27146469
#TreatL_S    -0.24592024 -0.07114642

# now compare adult total weights across treatments and generations in warmed populations only.
# this analysis compares fisheries selection treatments to fisheries size selection controls (uniform size selection) in warmed populations only

# load and prepare data for warmed 'F' generation analysis
Dev_weight_dataFH <- read.csv(file = "InputsZf/Weight_Dev_data_Warmed.csv")

Dev_weight_dataFH$Treat <- as.factor(paste(Dev_weight_dataFH$Temperature,Dev_weight_dataFH$Selection, sep = "_"))

Dev_weight_dataFH$Treat <- relevel(Dev_weight_dataFH$Treat , ref="H_C")

# first fit the full model to check fit
WeightFMH <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_weight_dataFH, REML = T)

plot(WeightFMH) # check resid spread 

qqnorm(residuals(WeightFMH))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=Dev_weight_dataFH$Generation,pearson=residuals(WeightFMH,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Treat=Dev_weight_dataFH$Treat,pearson=residuals(WeightFMH,type="pearson")),aes(x=Treat,y=pearson)) + geom_point() + theme_bw()

# now run the model selection for Deviance. 
# specify the full model
WeightFMHa <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_weight_dataFH, REML = F) 
# specify the additive model
WeightFMHb <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_weight_dataFH, REML = F) 

# test the models using anova
anova(WeightFMHa, WeightFMHb)
#Data: Dev_weight_dataFH
#Models:
# WeightFMHb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#WeightFMHa:  Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#.          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#WeightFMHb    8 682.37 724.03 -333.19   666.37                        
#WeightFMHa   10 673.90 725.98 -326.95   653.90 12.471  2   0.001959 **
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
WeightFMHbest <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_weight_dataFH, REML = T)

plot(allEffects(WeightFMHbest)) # gives good visual of parameter estimates

summary(WeightFMHbest)
#Fixed effects:
#                      Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept)         -7.048e-17  5.339e-02  8.274e+00   0.000   1.0000  
#TreatH_G             3.226e-02  4.862e-02  9.549e+00   0.664   0.5227  
#TreatH_S            -1.115e-01  4.862e-02  9.549e+00  -2.292   0.0460 *
#Generation          -7.890e-17  1.729e-02  5.093e+00   0.000   1.0000  
#TreatH_G:Generation  3.109e-02  1.445e-02  1.335e+03   2.152   0.0316 *
#TreatH_S:Generation -1.958e-02  1.445e-02  1.335e+03  -1.355   0.1756   

confint(WeightFMHbest)
#                           2.5 %       97.5 %
#.sig01               0.000000000  0.073772032
#.sig02               0.013731753  0.086670057
#.sig03               0.000000000  0.120293910
#.sigma               0.294976274  0.318211760
#(Intercept)         -0.097668496  0.097669193
#TreatH_G            -0.057538001  0.122066871
#TreatH_S            -0.201253557 -0.021648684
#Generation          -0.033572568  0.033572600
#TreatH_G:Generation  0.002772911  0.059400423
#TreatH_S:Generation -0.047890423  0.008737089


#### Juvenile length analysis ####

# Compare juvenile total lengths across treatments and generations. 

# load and prepare data for analysis during 'F' generations

JuvLen_data_dev <- read.csv(file = "InputsZf/JuvLength_Dev_data.csv")

Dev_Juvlength_dataF <- JuvLen_data_dev %>% filter (Generation > -1 & Generation < 6 )

Dev_Juvlength_dataF$Treat <- as.factor(paste(Dev_Juvlength_dataF$Temperature,Dev_Juvlength_dataF$Selection, sep = "_"))

Dev_Juvlength_dataF$Treat <- relevel(Dev_Juvlength_dataF$Treat , ref="L_C")

# first fit the full model to check fit
JuvLengthFM <- lmer(Deviance ~ Treat*Generation   + (1|Population) + (1|Tank) + (1|Generation), data = Dev_Juvlength_dataF, REML = T)# 

summary(JuvLengthFM)

plot(JuvLengthFM) # check resid spread 

qqnorm(residuals(JuvLengthFM))

# now run the model selection for Deviance. 
# specify the full model
JuvLengthFMa <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_Juvlength_dataF, REML = F) 
# specify the additive model
JuvLengthFMb <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_Juvlength_dataF, REML = F) 

# test the models using anova
anova(JuvLengthFMa, JuvLengthFMb)
#Data: Dev_Juvlength_dataF
#Models:
#JuvLengthFMb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#JuvLengthFMa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#             npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
#JuvLengthFMb   11 41554 41629 -20766    41532                         
#JuvLengthFMa   16 41384 41493 -20676    41352 179.55  5  < 2.2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
JuvLengthFMbest <- lmer(Deviance ~ Treat*Generation   + (1|Population) + (1|Tank) + (1|Generation), data = Dev_Juvlength_dataF, REML = T)# 

plot(allEffects(JuvLengthFMbest)) # gives good visual of parameter estimates

summary(JuvLengthFMbest)
#Fixed effects:
#                      Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)           -0.07365    0.85164    8.27778  -0.086   0.9331    
#TreatH_C               0.55182    0.70499   18.53643   0.783   0.4437    
#TreatH_G               0.05559    0.70499   18.53696   0.079   0.9380    
#TreatH_S               1.00055    0.70651   18.69434   1.416   0.1732    
#TreatL_G              -0.71572    0.69686   12.77393  -1.027   0.3234    
#TreatL_S               0.60885    0.69686   12.77393   0.874   0.3984    
#Generation             0.02009    0.24640    5.14675   0.082   0.9381    
#TreatH_C:Generation    0.85001    0.13193 6611.35375   6.443 1.25e-10 ***
#TreatH_G:Generation    1.22437    0.13193 6611.34755   9.280  < 2e-16 ***
#TreatH_S:Generation    0.18636    0.13726 6623.41670   1.358   0.1746    
#TreatL_G:Generation    0.61075    0.13250 6611.09422   4.609 4.11e-06 ***
#TreatL_S:Generation   -0.29009    0.13250 6611.09422  -2.189   0.0286 *    

confint(JuvLengthFMbest)
#                          2.5 %      97.5 %
#.sig01               0.32719258  0.88569495
#.sig02               0.00000000  0.75950492
#.sig03               0.46910731  1.60077841
#.sigma               5.33844201  5.52350934
#(Intercept)         -1.60733179  1.45841076
#TreatH_C            -0.69314512  1.79739200
#TreatH_G            -1.18999503  1.30056234
#TreatH_S            -0.24988199  2.24687788
#TreatL_G            -1.94061667  0.50886658
#TreatL_S            -0.61604246  1.83344079
#Generation          -0.44724784  0.48796200
#TreatH_C:Generation  0.59139688  1.10842517
#TreatH_G:Generation  0.96604046  1.48308120
#TreatH_S:Generation -0.07787766  0.46028358
#TreatL_G:Generation  0.35104370  0.87031690
#TreatL_S:Generation -0.54979978 -0.03052657

# CG1 analysis

# filter and prepare data for CG1 generation analysis
Dev_Juvlength_dataCG1 <- JuvLen_data_dev %>% filter (Generation == 6 )

Dev_Juvlength_dataCG1$Treat <- as.factor(paste(Dev_Juvlength_dataCG1$Temperature,Dev_Juvlength_dataCG1$Selection, sep = "_"))

Dev_Juvlength_dataCG1$Treat <- relevel(Dev_Juvlength_dataCG1$Treat , ref="L_C")

# first fit the full model to check fit
JuvLengthCG1M <- lmer(Deviance ~ Treat  + (1|Tank), data = Dev_Juvlength_dataCG1, REML = T)

summary(JuvLengthCG1M)

plot(JuvLengthCG1M) # check resid spread 

qqnorm(residuals(JuvLengthCG1M))

# now run the model selection for Deviance. 
# specify the full model
JuvLengthCG1Ma <- lmer(Deviance ~  Treat  + (1|Tank), data = Dev_Juvlength_dataCG1, REML = F) 
# specify the null
JuvLengthCG1Mb <- lmer(Deviance ~   + (1|Tank), data = Dev_Juvlength_dataCG1, REML = F) 

# test the models using anova
anova(JuvLengthCG1Ma, JuvLengthCG1Mb)
#Data: Dev_Juvlength_dataCG1
#Models:
#JuvLengthCG1Mb: Deviance ~ +(1 | Tank)
#JuvLengthCG1Ma: Deviance ~ Treat + (1 | Tank)
#               npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
#JuvLengthCG1Mb    3 6580.9 6595.7 -3287.4   6574.9                       
#JuvLengthCG1Ma    8 6576.6 6616.1 -3280.3   6560.6 14.258  5    0.01405 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
JuvLengthCG1Mbest <- lmer(Deviance ~  Treat  + (1|Tank) , data = Dev_Juvlength_dataCG1, REML = T) 

plot(allEffects(JuvLengthCG1Mbest)) # gives good visual of parameter estimates

summary(JuvLengthCG1Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept) -2.547e-13  4.509e-01  1.014e+03   0.000   1.0000  
#TreatH_C    -9.111e-01  6.377e-01  1.014e+03  -1.429   0.1534  
#TreatH_G     2.056e-01  6.377e-01  1.014e+03   0.322   0.7473  
#TreatH_S     7.944e-01  7.130e-01  1.014e+03   1.114   0.2654  
#TreatL_G     1.228e+00  6.377e-01  1.014e+03   1.925   0.0545 .
#TreatL_S     8.667e-01  6.377e-01  1.014e+03   1.359   0.1744     

confint(JuvLengthCG1Mbest)
#                  2.5 %    97.5 %
#.sig01       0.00000000 0.4888603
#.sigma       5.77931753 6.3034201
#(Intercept) -0.88195653 0.8819566
#TreatH_C    -2.15838600 0.3361638
#TreatH_G    -1.04171933 1.4528305
#TreatH_S    -0.60013355 2.1892132
#TreatL_G    -0.01957581 2.4751314
#TreatL_S    -0.38068692 2.1140203

# CG2 analysis

# filter and prepare data for CG2 generation analysis
Dev_Juvlength_dataCG2 <- JuvLen_data_dev %>% filter (Generation > 6 )

Dev_Juvlength_dataCG2$Treat <- as.factor(paste(Dev_Juvlength_dataCG2$Temperature,Dev_Juvlength_dataCG2$Selection, sep = "_"))

Dev_Juvlength_dataCG2$Treat <- relevel(Dev_Juvlength_dataCG2$Treat , ref="L_C")

# first fit the full model to check fit
JuvLengthCG2M <- lmer(Deviance ~ Treat + (1|Tank), data = Dev_Juvlength_dataCG2, REML = T)

summary(JuvLengthCG2M)

plot(JuvLengthCG2M) # check resid spread 

qqnorm(residuals(JuvLengthCG2M))

# now run the model selection for Deviance. 
# specify the full model
JuvLengthCG2Ma <- lmer(Deviance ~ Treat  + (1|Tank), data = Dev_Juvlength_dataCG2, REML = F) 
# specify the null model
JuvLengthCG2Mb <- lmer(Deviance ~   +  (1|Tank), data = Dev_Juvlength_dataCG2, REML = F) 

# test the models using anova
anova(JuvLengthCG2Ma, JuvLengthCG2Mb)
#Data: Dev_Juvlength_dataCG2
#Models:
#JuvLengthCG2Mb: Deviance ~ +(1 | Tank)
#JuvLengthCG2Ma: Deviance ~ Treat + (1 | Tank)
#               npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#JuvLengthCG2Mb    3 6255.1 6269.9 -3124.6   6249.1                     
#JuvLengthCG2Ma    8 6260.5 6299.9 -3122.2   6244.5 4.6541  5     0.4595

# now fit the supported model
JuvLengthCG2Mbest <- lmer(Deviance ~  + (1|Tank) , data = Dev_Juvlength_dataCG2, REML = T) 

#plot(allEffects(JuvLengthCG2Mbest)) # model containts no terms other than the constant

summary(JuvLengthCG2Mbest)
#Fixed effects:
#             Estimate Std. Error        df t value Pr(>|t|)
#(Intercept)    0.1614     0.1622 1019.0000   0.995     0.32  

confint(JuvLengthCG2Mbest)
#                 2.5 %    97.5 %
#.sig01       0.0000000 0.5447584
#.sigma       4.9608916 5.4107757
#(Intercept) -0.1566111 0.4794870

# now compare juvenile total lengths across treatments and generations in warmed populations only.
# this analysis compares fisheries selection treatments to fisheries size selection controls (uniform size selection) in warmed populations only

# load and prepare data for warmed 'F' generation analysis
Dev_JuvLength_dataFH <- read.csv(file = "InputsZf/JuvLength_Dev_data_Warmed.csv")

Dev_JuvLength_dataFH$Treat <- as.factor(paste(Dev_JuvLength_dataFH$Temperature,Dev_JuvLength_dataFH$Selection, sep = "_"))

Dev_JuvLength_dataFH$Treat <- relevel(Dev_JuvLength_dataFH$Treat , ref="H_C")

# first fit the full model to check fit
JuvLengthFMH <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_JuvLength_dataFH, REML = T)# 

plot(JuvLengthFMH) # check resid spread 

qqnorm(residuals(JuvLengthFMH))

# now run the model selection for Deviance. 
# specify the full model
JuvLengthFMHa <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_JuvLength_dataFH, REML = F) 
# specify the additive model
JuvLengthFMHb <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_JuvLength_dataFH, REML = F) 

# test the models using anova
anova(JuvLengthFMHa, JuvLengthFMHb)
#Data: Dev_JuvLength_dataFH
#Models:
#JuvLengthFMHb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#JuvLengthFMHa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#              npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
#JuvLengthFMHb    8 20531 20580 -10258    20515                         
#JuvLengthFMHa   10 20477 20538 -10228    20457 58.473  2  2.008e-13 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
JuvLengthFMHbest <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_JuvLength_dataFH, REML = T)# 

plot(allEffects(JuvLengthFMHbest)) # gives good visual of parameter estimates

summary(JuvLengthFMHbest)
#Fixed effects:
#                     Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)            0.07111    0.85415    8.60923   0.083  0.93556    
#TreatH_G              -0.49623    0.74047    5.85704  -0.670  0.52828    
#TreatH_S               0.44744    0.74193    5.90297   0.603  0.56889    
#Generation            -0.01939    0.23276    4.95740  -0.083  0.93686    
#TreatH_G:Generation    0.37415    0.12996 3283.14028   2.879  0.00402 ** 
#TreatH_S:Generation   -0.66285    0.13542 3290.26752  -4.895 1.03e-06 ***

confint(JuvLengthFMHbest)
#                         2.5 %     97.5 %
#.sig01               0.3050417  1.2362371
#.sig02               0.4242169  1.5687936
#.sig03               0.0000000  1.5284403
#.sigma               5.2221605  5.4809441
#(Intercept)         -1.5009123  1.6440245
#TreatH_G            -1.8742877  0.8823196
#TreatH_S            -0.9334939  1.8279828
#Generation          -0.4737346  0.4345062
#TreatH_G:Generation  0.1196674  0.6291002
#TreatH_S:Generation -0.9247433 -0.3936081

#### Age at maturity (A50) analysis ####

# Compare age at maturity (A50) across treatments and generations. 

# load and prepare data for analysis during 'F' generations
A50_data_dev <- read.csv(file = "InputsZf/A50_Dev_data.csv")

Dev_A50_Fdata <- A50_data_dev %>% filter (Generation > -1 & Generation < 6 )

Dev_A50_Fdata$Treat <- as.factor(paste(Dev_A50_Fdata$Temperature,Dev_A50_Fdata$Selection, sep = "_"))

Dev_A50_Fdata$Treat <- relevel(Dev_A50_Fdata$Treat , ref="L_C")

# first fit the full model to check fit
A50FM <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_Fdata, REML = T)# 

plot(A50FM) # check resid spread 

qqnorm(residuals(A50FM))

# now run the model selection for Deviance. 
# specify the full model
A50FM1a <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_Fdata, REML = F) 
# specify the model with 2-way interactions
A50FM1b <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_Fdata, REML = F) 

# test the models using anova
anova(A50FM1a, A50FM1b)
#Data: Dev_A50_Fdata
#Models:
#A50FM1b: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#A50FM1a: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#        npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#A50FM1b   11 663.71 693.01 -320.85   641.71                         
#A50FM1a   16 648.99 691.61 -308.50   616.99 24.716  5  0.0001581 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
A50FM1best <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_Fdata, REML = T)# 

plot(allEffects(A50FM1best)) # gives good visual of parameter estimates

summary(A50FM1best)
#Fixed effects:
#                      Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept)          7.826e-14  2.057e+00  4.357e+01   0.000   1.0000  
#TreatH_C            -1.871e+00  2.827e+00  5.839e+01  -0.662   0.5107  
#TreatH_G            -3.554e+00  2.827e+00  5.839e+01  -1.257   0.2136  
#TreatH_S            -1.349e+00  2.846e+00  5.899e+01  -0.474   0.6372  
#TreatL_G            -1.658e+00  2.967e+00  6.265e+01  -0.559   0.5784  
#TreatL_S            -1.434e+00  2.827e+00  5.839e+01  -0.507   0.6140  
#Generation          -2.081e-14  6.569e-01  4.875e+01   0.000   1.0000  
#TreatH_C:Generation -2.286e+00  9.008e-01  7.813e+01  -2.538   0.0131 *
#TreatH_G:Generation -1.831e+00  9.008e-01  7.813e+01  -2.032   0.0455 *
#TreatH_S:Generation -1.189e+00  9.338e-01  7.962e+01  -1.274   0.2065  
#TreatL_G:Generation  1.186e+00  9.338e-01  7.962e+01   1.270   0.2076  
#TreatL_S:Generation  1.036e+00  9.008e-01  7.813e+01   1.150   0.2537   

confint(A50FM1best)
#                         2.5 %     97.5 %
#.sig01               0.0000000  1.2764584
#.sig02               0.0000000  1.8605133
#.sigma               3.8907773  5.1157944
#(Intercept)         -3.6821790  3.6821892
#TreatH_C            -7.0636651  3.3223507
#TreatH_G            -8.7472852  1.6387305
#TreatH_S            -6.5591383  3.8992658
#TreatL_G            -7.1059882  3.8169652
#TreatL_S            -6.6267973  3.7595929
#Generation          -1.2165162  1.2165175
#TreatH_C:Generation -4.0013852 -0.5708734
#TreatH_G:Generation -3.5459644 -0.1154526
#TreatH_S:Generation -2.9804981  0.5722804
#TreatL_G:Generation -0.5934836  2.9586987
#TreatL_S:Generation -0.6795324  2.7509793

# CG1 analysis

# filter and prepare data for CG1 generation analysis
A50_Dev_dataCG1 <- A50_data_dev %>% filter (Generation == 6 )

A50_Dev_dataCG1$Treat <- as.factor(paste(A50_Dev_dataCG1$Temperature,A50_Dev_dataCG1$Selection, sep = "_"))

A50_Dev_dataCG1$Treat <- relevel(A50_Dev_dataCG1$Treat , ref="L_C")

# first fit the full model to check fit
A50CG1M <- lmer(Deviance ~ Treat + (1|Tank), data = A50_Dev_dataCG1, REML = T)

summary(A50CG1M)

plot(A50CG1M) # check resid spread 

qqnorm(residuals(A50CG1M))

# now run the model selection for Deviance. 
# specify the full model
A50CG1Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = A50_Dev_dataCG1, REML = F) 
# specify the null model
A50CG1Mb <- lmer(Deviance ~   + (1|Tank), data = A50_Dev_dataCG1, REML = F) 

# test the models using anova
anova(A50CG1Ma, A50CG1Mb)
#Data: A50_Dev_dataCG1
#Models:
#A50CG1Mb: Deviance ~ +(1 | Tank)
#A50CG1Ma: Deviance ~ Treat + (1 | Tank)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
#A50CG1Mb    3 110.53 113.03 -52.267  104.533                       
#A50CG1Ma    8 109.55 116.21 -46.774   93.547 10.986  5    0.05166 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
A50CG1Mbest <- lmer(Deviance ~   + (1|Tank) , data = A50_Dev_dataCG1, REML = T) 

#plot(allEffects(A50CG1Mbest)) # the model contains no fixed terms so don't plot

summary(A50CG1Mbest)
#Fixed effects:
#            Estimate Std. Error     df t value Pr(>|t|)
#(Intercept)   -1.110      1.309 16.000  -0.848    0.409

confint(A50CG1Mbest)
#                2.5 %   97.5 %
#.sig01       0.000000 4.082721
#.sigma       3.869778 7.644369
#(Intercept) -3.745531 1.539833

# CG2 analysis

# filter and prepare data for CG2 generation analysis
A50_Dev_dataCG2 <- A50_data_dev %>% filter (Generation > 6 )

A50_Dev_dataCG2$Treat <- as.factor(paste(A50_Dev_dataCG2$Temperature,A50_Dev_dataCG2$Selection, sep = "_"))

A50_Dev_dataCG2$Treat <- relevel(A50_Dev_dataCG2$Treat , ref="L_C")

# first fit the full model to check fit
A50CG2M <- lmer(Deviance ~ Treat  + (1|Tank), data = A50_Dev_dataCG2, REML = T)

summary(A50CG2M)

plot(A50CG2M) # check resid spread 

qqnorm(residuals(A50CG2M))

# check which residuals are greater than 5 or smaller than -5 
which(residuals(A50CG2M) > 5)
which(residuals(A50CG2M) < -5)
# remove them as they look to be outliers
A50_Dev_dataCG2a <- A50_Dev_dataCG2[-which(residuals(A50CG2M) > 5),]
A50_Dev_dataCG2a <- A50_Dev_dataCG2a[-which(residuals(A50CG2M) < -5),]

# fit the model again and check fit
A50CG2M1 <- lmer(Deviance ~ Treat + (1|Tank), data = A50_Dev_dataCG2a, REML = T)

plot(A50CG2M1) # check resid spread 

qqnorm(residuals(A50CG2M1))

# now run the model selection for Deviance. 
# specify the full model
A50CG2M1a <- lmer(Deviance ~  Treat  + (1|Tank), data = A50_Dev_dataCG2a, REML = F) 
# specify the model with 2-way interactions
A50CG2M1b <- lmer(Deviance ~   + (1|Tank), data = A50_Dev_dataCG2a, REML = F) 

# test the models using anova
anova(A50CG2M1a, A50CG2M1b)
#Data: A50_Dev_dataCG2a
#Models:
#A50CG2M1b: Deviance ~ +(1 | Tank)
#A50CG2M1a: Deviance ~ Treat + (1 | Tank)
#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#A50CG2M1b    3 72.895 75.019 -33.447   66.895                     
#A50CG2M1a    8 75.487 81.152 -29.744   59.487 7.4074  5     0.1921

# now fit the supported model
A50_CG2Mbest <- lmer(Deviance ~ + (1|Tank) , data = A50_Dev_dataCG2a, REML = T) 

#plot(allEffects(A50_CG2Mbest)) # model contains no terms

summary(A50_CG2Mbest)
#Fixed effects:
#            Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)  -0.2767     0.6013 14.0000   -0.46    0.652

confint(A50_CG2Mbest)
#
#                2.5 %    97.5 %
#.sig01       0.000000 2.1945490
#.sigma       1.634104 3.3770357
#(Intercept) -1.492274 0.9387795

# now compare age at maturity (A50) across treatments and generations in warmed populations only.
# this analysis compares fisheries selection treatments to fisheries size selection controls (uniform size selection) in warmed populations only

# load and prepare data for warmed 'F' generation analysis
Dev_A50_dataFH <- read.csv(file = "InputsZf/A50_Dev_data_Warmed.csv")

Dev_A50_dataFH$Treat <- as.factor(paste(Dev_A50_dataFH$Temperature,Dev_A50_dataFH$Selection, sep = "_"))

Dev_A50_dataFH$Treat <- relevel(Dev_A50_dataFH$Treat , ref="H_C")

# first fit the full model and check fit
A50FMH <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_dataFH, REML = T)# 

summary(A50FMH)

plot(A50FMH) # check resid spread 

qqnorm(residuals(A50FMH))

# now run the model selection for Deviance. 
# specify the full model
A50FMHa <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_dataFH, REML = F) 
# specify the additive model
A50FMHb <- lmer(Deviance ~ Treat + Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_dataFH, REML = F) 

# test the models using anova
anova(A50FMHa, A50FMHb)
#Data: Dev_A50_dataFH
#Models:
#A50FMHb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#A50FMHa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#        npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#A50FMHb    8 324.56 340.32 -154.28   308.56                     
#A50FMHa   10 327.54 347.24 -153.77   307.54 1.0171  2     0.6014

# now fit the supported model
A50FMHbest <- lmer(Deviance ~ Treat + Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_A50_dataFH, REML = T)# 

plot(allEffects(A50FMHbest)) # gives good visual of parameter estimates

summary(A50FMHbest)
#Fixed effects:
#                 Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)       -1.1281     1.7186  6.2895  -0.656    0.535
#TreatH_G          -0.5451     1.6715  5.3312  -0.326    0.757
#TreatH_S           3.0319     1.6926  5.5806   1.791    0.127
#Generation         0.4512     0.4599  4.0463   0.981    0.382

confint(A50FMHbest)
#                  2.5 %   97.5 %
#.sig01       0.00000000 2.319769
#.sig02       0.00000000 2.720996
#.sig03       0.00000000 2.028918
#.sigma       3.63745242 5.446785
#(Intercept) -4.10733432 1.894314
#TreatH_G    -3.57658089 2.486452
#TreatH_S    -0.08260761 6.069673
#Generation  -0.41768855 1.305441

#### Length at maturity (L50) analysis ####

# Compare length at maturity (L50) across treatments and generations. 

# load and prepare data for analysis during 'F' generations
L50_data_dev <- read.csv(file = "InputsZf/L50_Dev_data.csv")

Dev_L50_Fdata <- L50_data_dev %>% filter (Generation > -1 & Generation < 6 )

Dev_L50_Fdata$Treat <- as.factor(paste(Dev_L50_Fdata$Temperature,Dev_L50_Fdata$Selection, sep = "_"))

Dev_L50_Fdata$Treat <- relevel(Dev_L50_Fdata$Treat , ref="L_C")

# first fit the full model to check fit
L50FM <- lmer(Deviance ~ Treat*Generation + (1|Population)  + (1|Tank) + (1|Generation), data = Dev_L50_Fdata, REML = T)

plot(L50FM) # check resid spread 

qqnorm(residuals(L50FM))

# now run the model selection for Deviance. 
# specify the full model
L50FMa <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_Fdata, REML = F) 
# specify the additive model
L50FMb <- lmer(Deviance ~ Treat + Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_Fdata, REML = F) 

# test the models using anova
anova(L50FMa, L50FMb)
#Data: Dev_L50_Fdata
#Models:
#L50FMb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#L50FMa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#       npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#L50FMb   11 467.39 496.80 -222.70   445.39                        
#L50FMa   16 460.06 502.82 -214.03   428.06 17.338  5   0.003901 **
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
L50FMbest <- lmer(Deviance ~ Treat*Generation + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_Fdata, REML = T)# 

plot(allEffects(L50FMbest)) # gives good visual of parameter estimates

summary(L50FMbest)
#                      Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept)          9.044e-16  9.580e-01  2.190e+01   0.000   1.0000  
#TreatH_C            -4.633e-01  1.150e+00  4.839e+01  -0.403   0.6887  
#TreatH_G            -1.755e+00  1.150e+00  4.839e+01  -1.527   0.1333  
#TreatH_S            -1.755e-01  1.157e+00  4.903e+01  -0.152   0.8800  
#TreatL_G            -6.699e-01  1.116e+00  3.943e+01  -0.600   0.5519  
#TreatL_S            -1.231e+00  1.116e+00  3.943e+01  -1.102   0.2770  
#Generation           6.513e-16  2.955e-01  1.953e+01   0.000   1.0000  
#TreatH_C:Generation -3.214e-01  3.443e-01  7.916e+01  -0.933   0.3534  
#TreatH_G:Generation  2.443e-01  3.443e-01  7.916e+01   0.709   0.4802  
#TreatH_S:Generation -5.293e-01  3.575e-01  8.027e+01  -1.481   0.1426  
#TreatL_G:Generation  7.890e-01  3.443e-01  7.916e+01   2.291   0.0246 *
#TreatL_S:Generation  2.896e-01  3.443e-01  7.916e+01   0.841   0.4028  


confint(L50FMbest)
#                          2.5 %    97.5 %
#.sig01               0.0000000 0.9747317
#.sig02               0.0868717 1.2983474
#.sigma               1.4990427 1.9924417
#(Intercept)         -1.6922631 1.6922619
#TreatH_C            -2.5419535 1.6152899
#TreatH_G            -3.8339199 0.3233235
#TreatH_S            -2.2578060 1.9275588
#TreatL_G            -2.6781341 1.3384219
#TreatL_S            -3.2388961 0.7776600
#Generation          -0.5450036 0.5450035
#TreatH_C:Generation -0.9847221 0.3419027
#TreatH_G:Generation -0.4190483 0.9075765
#TreatH_S:Generation -1.2234262 0.1523546
#TreatL_G:Generation  0.1256637 1.4522885
#TreatL_S:Generation -0.3737019 0.9529230

# CG1 analysis

# filter and prepare data for CG1 generation analysis
L50_Dev_dataCG1 <- L50_data_dev %>% filter (Generation == 6 )

L50_Dev_dataCG1$Treat <- as.factor(paste(L50_Dev_dataCG1$Temperature,L50_Dev_dataCG1$Selection, sep = "_"))

L50_Dev_dataCG1$Treat <- relevel(L50_Dev_dataCG1$Treat , ref="L_C")

# first fit the full model to check fit
L50CG1M <- lmer(Deviance ~ Treat + (1|Tank), data = L50_Dev_dataCG1, REML = T)

summary(L50CG1M)

plot(L50CG1M) # check resid spread 

qqnorm(residuals(L50CG1M))

# now run the model selection for Deviance. 
# specify the full model
L50CG1Ma <- lmer(Deviance ~  Treat + (1|Tank), data = L50_Dev_dataCG1, REML = F) 
# specify the null model
L50CG1Mb <- lmer(Deviance ~   + (1|Tank), data = L50_Dev_dataCG1, REML = F) 

# test the models using anova
anova(L50CG1Ma, L50CG1Mb)
#Data: L50_Dev_dataCG1
#Models:
#L50CG1Mb: Deviance ~ +(1 | Tank)
#L50CG1Ma: Deviance ~ Treat + (1 | Tank)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#L50CG1Mb    3 70.940 73.440 -32.470   64.940                        
#L50CG1Ma    8 65.027 71.693 -24.514   49.027 15.912  5   0.007098 **
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
L50CG1Mbest <- lmer(Deviance ~  Treat  + (1|Tank) , data = L50_Dev_dataCG1, REML = T) 

plot(allEffects(L50CG1Mbest)) # gives good visual of parameter estimates

summary(L50CG1Mbest)
#Fixed effects:
#              Estimate Std. Error         df t value Pr(>|t|)  
#(Intercept)  1.077e-16  7.345e-01  1.100e+01   0.000   1.0000  
#TreatH_C     3.786e-01  1.039e+00  1.100e+01   0.364   0.7224  
#TreatH_G     1.843e+00  1.039e+00  1.100e+01   1.774   0.1037  
#TreatH_S    -2.548e+00  1.161e+00  1.100e+01  -2.194   0.0506 .
#TreatL_G     1.051e+00  1.039e+00  1.100e+01   1.012   0.3332  
#TreatL_S     1.295e+00  1.039e+00  1.100e+01   1.247   0.2384   

confint(L50CG1Mbest)
#                 2.5 %     97.5 %
#.sig01       0.0000000  0.7295195
#.sigma       0.7563232  1.4939839
#(Intercept) -1.2265471  1.2265473
#TreatH_C    -1.3560415  2.1131579
#TreatH_G     0.1080417  3.5772411
#TreatH_S    -4.4873359 -0.6086240
#TreatL_G    -0.6831736  2.7860656
#TreatL_S    -0.4396732  3.0295661

# CG2 analysis 

# filter and prepare data for CG2 generation analysis
L50_Dev_dataCG2 <- L50_data_dev %>% filter (Generation > 6 )

L50_Dev_dataCG2$Treat <- as.factor(paste(L50_Dev_dataCG2$Temperature,L50_Dev_dataCG2$Selection, sep = "_"))

L50_Dev_dataCG2$Treat <- relevel(L50_Dev_dataCG2$Treat , ref="L_C")

# first fit the full model and check fit
L50CG2M <- lmer(Deviance ~ Treat  + (1|Tank), data = L50_Dev_dataCG2, REML = T)

summary(L50CG2M)

plot(L50CG2M) # check resid spread 

qqnorm(residuals(L50CG2M))

# now run the model selection for Deviance. 
# specify the full model
L50CG2Ma <- lmer(Deviance ~  Treat  + (1|Tank), data = L50_Dev_dataCG2, REML = F) 
# specify the null model
L50CG2Mb <- lmer(Deviance ~   + (1|Tank), data = L50_Dev_dataCG2, REML = F) 

# test the models using anova
anova(L50CG2Ma, L50CG2Mb)
#Data: L50_Dev_dataCG2
#Models:
#L50CG2Mb: Deviance ~ +(1 | Tank)
#L50CG2Ma: Deviance ~ Treat + (1 | Tank)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#L50CG2Mb    3 76.734 79.234 -35.367   70.734                     
#L50CG2Ma    8 77.929 84.595 -30.965   61.929 8.8054  5     0.1171

# now fit the supported model
L50CG2Mbest <- lmer(Deviance ~   + (1|Tank) , data = L50_Dev_dataCG2, REML = T) 

#plot(allEffects(L50CG2Mbest)) # the model containts no fixed terms so don't plot

summary(L50CG2Mbest)
#Fixed effects:
#            Estimate Std. Error     df t value Pr(>|t|)
#(Intercept)   0.8694     0.6260 5.2301   1.389    0.221

confint(L50CG2Mbest)
#                2.5 %   97.5 %
#.sig01       0.000000 2.588093
#.sigma       1.211645 2.761081
#(Intercept) -0.451743 2.194177

# now compare length at maturity (L50) across treatments and generations in warmed populations only.
# this analysis compares fisheries selection treatments to fisheries size selection controls (uniform size selection) in warmed populations only

# load and prepare data for warmed 'F' generation analysis
Dev_L50_dataFH <- read.csv(file = "InputsZf/L50_Dev_data_Warmed.csv")

Dev_L50_dataFH$Treat <- as.factor(paste(Dev_L50_dataFH$Temperature,Dev_L50_dataFH$Selection, sep = "_"))

Dev_L50_dataFH$Treat <- relevel(Dev_L50_dataFH$Treat , ref="H_C")

# first fit the full model and check fit
L50FMH <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_dataFH, REML = T)# 

plot(L50FMH) # check resid spread 

qqnorm(residuals(L50FMH))

# now run the model selection for Deviance. 
# specify the full model
L50FMHa <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_dataFH, REML = F) 
# specify the additive model
L50FMHb <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_dataFH, REML = F) 

# test the models using anova
anova(L50FMHa, L50FMHb)
#Data: Dev_L50_dataFH
#Models:
#L50FMHb: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#L50FMHa: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#        npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
#L50FMHb    8 232.35 248.11 -108.17   216.35                       
#L50FMHa   10 231.19 250.89 -105.59   211.19 5.1622  2    0.07569 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
L50FMHbest <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_L50_dataFH, REML = T)# 

plot(allEffects(L50FMHbest)) # gives good visual of parameter estimates

summary(L50FMHbest)
#Fixed effects:
#                 Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)  -0.3116     0.6085 10.0911  -0.512    0.620
#TreatH_G      0.1222     0.6449 44.9776   0.190    0.851
#TreatH_S     -0.2131     0.6548 45.1318  -0.325    0.746
#Generation    0.1247     0.1600  4.1817   0.779    0.478

confint(L50FMHbest)
#                 2.5 %    97.5 %
#.sig01       0.0000000 0.7745499
#.sig02       0.0000000 0.9240166
#.sig03       0.0000000 0.8285140
#.sigma       1.5574889 2.2826698
#(Intercept) -1.4686000 0.8470772
#TreatH_G    -1.1158748 1.3603106
#TreatH_S    -1.4711077 1.0429200
#Generation  -0.1783164 0.4270886


#### Adult female GSI analysis ####

# Compare adult female gonado-somatic index (GSI) across treatments and generations. 

# load and prepare data for analysis during 'F' generations
GSI_data_dev <- read.csv(file = "InputsZf/GSI_Dev_data.csv")

GSI_dev_dataF <- GSI_data_dev %>% filter (Generation > -1 & Generation < 5 )

GSI_dev_dataF$Treat <- as.factor(paste(GSI_dev_dataF$Temperature,GSI_dev_dataF$Selection, sep = "_"))

GSI_dev_dataF$Treat <- relevel(GSI_dev_dataF$Treat , ref="L_C")

# first fit the full model and check fit
GSIFM <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = GSI_dev_dataF, REML = T)

plot(GSIFM) # check resid spread 

qqnorm(residuals(GSIFM))

# now run the model selection for Deviance. 
# specify the full model
GSIFM1a <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = GSI_dev_dataF, REML = F) 
# specify the additive model
GSIFM1b <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = GSI_dev_dataF, REML = F) 

# test the models using anova
anova(GSIFM1a, GSIFM1b)
#Data: GSI_dev_dataF1
#Models:
#GSIFM1b: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#GSIFM1a: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#        npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#GSIFM1b   11 4380.8 4430.3 -2179.4   4358.8                     
#GSIFM1a   16 4384.8 4456.8 -2176.4   4352.8 6.0299  5     0.3033

# now fit the supported model
GSIFM1best <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = GSI_dev_dataF, REML = T)# 

plot(allEffects(GSIFM1best)) # gives good visual of parameter estimates

summary(GSIFM1best)
#Fixed effects:
#            Estimate Std. Error       df t value Pr(>|t|)   
#(Intercept)   2.9825     1.2480   4.1545   2.390  0.07279 . 
#TreatH_C     -1.3475     0.8551 657.2009  -1.576  0.11553   
#TreatH_G     -2.7649     0.8439 657.1886  -3.276  0.00111 **
#TreatH_S     -1.7127     0.8945 657.2809  -1.915  0.05596 . 
#TreatL_G     -1.1763     0.8806 657.6846  -1.336  0.18208   
#TreatL_S     -1.0851     0.9308 657.7980  -1.166  0.24409   
#Generation   -1.0577     0.3897   2.2249  -2.714  0.10097     
   

confint(GSIFM1best)
#                 2.5 %      97.5 %
#.sig01       0.0000000  0.70583304
#.sig02       0.0000000  0.80102083
#.sig03       0.0000000  1.42501337
#.sigma       6.0440452  6.73136731
#(Intercept)  0.7706504  5.19841732
#TreatH_C    -3.0259329  0.31749946
#TreatH_G    -4.4076879 -1.10801569
#TreatH_S    -3.4508975  0.04685041
#TreatL_G    -2.9096011  0.53342926
#TreatL_S    -2.9260745  0.71403399
#Generation  -1.7678767 -0.34789583

# Lets look at F5 generation to see the effects after treatments were relaxed

# filter and prepare data for F5 generation analysis
GSI_data_devF5 <- GSI_data_dev %>% filter (Generation == 5 )

GSI_data_devF5$Treat <- as.factor(paste(GSI_data_devF5$Temperature,GSI_data_devF5$Selection, sep = "_"))

GSI_data_devF5$Treat <- relevel(GSI_data_devF5$Treat , ref="L_C")

# first fit the full model to check fit
GSIF5M <- lmer(Deviance ~ Treat  + (1|Tank), data = GSI_data_devF5, REML = T)

summary(GSIF5M)

plot(GSIF5M) # check resid spread 

qqnorm(residuals(GSIF5M))

# now run the model selection for Deviance. 
# specify the full model
GSIF5Ma <- lmer(Deviance ~  Treat   + (1|Tank), data = GSI_data_devF5, REML = F) 
# specify the null model
GSIF5Mb <- lmer(Deviance ~     + (1|Tank), data = GSI_data_devF5, REML = F) 

# test the models using anova
anova(GSIF5Ma, GSIF5Mb)
#Data: GSI_data_devF5
#Models:
#GSIF5Mb: Deviance ~ +(1 | Tank)
#GSIF5Ma: Deviance ~ Treat + (1 | Tank)
#        npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#GSIF5Mb    3 1280.4 1290.2 -637.22   1274.4                     
#GSIF5Ma    8 1285.2 1311.3 -634.58   1269.2 5.2858  5      0.382

# now fit the supported model
GSI_F5Mbest <- lmer(Deviance ~  + (1|Tank) , data = GSI_data_devF5, REML = T) 

#plot(allEffects(GSI_F5Mbest)) # model containts no fixed terms to don't plot

summary(GSI_F5Mbest)
#Fixed effects:
#            Estimate Std. Error      df t value Pr(>|t|)  
#(Intercept)  -1.8316     0.6311  3.9800  -2.902   0.0443 *

confint(GSI_F5Mbest)
#                2.5 %     97.5 %
#.sig01       0.000000  2.6688396
#.sigma       5.822464  7.1351032
#(Intercept) -3.172513 -0.4754628

# CG1 analsis 

# filter and prepare data for CG1 generation analysis
GSI_data_devCG1 <- GSI_data_dev %>% filter (Generation == 6 )

GSI_data_devCG1$Treat <- as.factor(paste(GSI_data_devCG1$Temperature,GSI_data_devCG1$Selection, sep = "_"))

GSI_data_devCG1$Treat <- relevel(GSI_data_devCG1$Treat , ref="L_C")

# first fit the full model to check fit
GSICG1M <- lmer(Deviance ~ Treat  + (1|Tank), data = GSI_data_devCG1, REML = T)

summary(GSICG1M)

plot(GSICG1M) # check resid spread 

qqnorm(residuals(GSICG1M))

# now run the model selection for Deviance. 
# specify the full model
GSICG1Ma <- lmer(Deviance ~  Treat  + (1|Tank), data = GSI_data_devCG1, REML = F) 
# specify the null model
GSICG1Mb <- lmer(Deviance ~    + (1|Tank), data = GSI_data_devCG1, REML = F) 

# test the models using anova
anova(GSICG1Ma, GSICG1Mb)
#Data: GSI_data_devCG1
#Models:
#GSICG1Mb: Deviance ~ +(1 | Tank)
#GSICG1Ma: Deviance ~ Treat + (1 | Tank)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#GSICG1Mb    3 1423.7 1433.7 -708.87   1417.7                     
#GSICG1Ma    8 1431.3 1457.9 -707.64   1415.3 2.4699  5      0.781

# now fit the supported model
GSICG1Mbest <- lmer(Deviance ~  + (1|Tank) , data = GSI_data_devCG1, REML = T) 

#plot(allEffects(GSICG1Mbest)) # model containts no fixed terms to don't plot

summary(GSICG1Mbest)
#Fixed effects:
#            Estimate Std. Error       df t value Pr(>|t|)
#(Intercept)   0.6982     0.5177 206.0000   1.349    0.179

confint(GSICG1Mbest)
#                2.5 %   97.5 %
#.sig01       0.0000000 1.668519
#.sigma       6.7681317 8.207743
#(Intercept) -0.3187079 1.715116

# CG2 analysis

# filter and prepare data for CG2 generation analysis
GSI_data_devCG2 <- GSI_data_dev %>% filter (Generation == 7 )

GSI_data_devCG2$Treat <- as.factor(paste(GSI_data_devCG2$Temperature,GSI_data_devCG2$Selection, sep = "_"))

GSI_data_devCG2$Treat <- relevel(GSI_data_devCG2$Treat , ref="L_C")

# first fit the full model to check fit
GSICG2M <- lmer(Deviance ~ Treat  + (1|Tank), data = GSI_data_devCG2, REML = T)

summary(GSICG2M)

plot(GSICG2M) # check resid spread 

qqnorm(residuals(GSICG2M))

# now run the model selection for Deviance. 
# specify the full model
GSICG2Ma <- lmer(Deviance ~  Treat + (1|Tank), data = GSI_data_devCG2, REML = F) 
# specify the null model
GSICG2Mb <- lmer(Deviance ~     + (1|Tank), data = GSI_data_devCG2, REML = F) 

# test the models using anova
anova(GSICG2Ma, GSICG2Mb)
#Data: GSI_data_devCG2
#Models:
#GSICG2Mb: Deviance ~ +(1 | Tank)
#GSICG2Ma: Deviance ~ Treat + (1 | Tank)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#GSICG2Mb    3 1232.3 1241.9 -613.14   1226.3                     
#GSICG2Ma    8 1236.2 1262.0 -610.12   1220.2 6.0471  5     0.3017

# now fit the supported model
GSICG2Mbest <- lmer(Deviance ~ + (1|Tank) , data = GSI_data_devCG2, REML = T) 

#plot(allEffects(GSICG2Mbest)) # model containts no fixed terms to don't plot

summary(GSICG2Mbest)
#Fixed effects:
#            Estimate Std. Error       df t value Pr(>|t|)  
#(Intercept)  -1.2468     0.5009 183.0000  -2.489   0.0137 *

confint(GSICG2Mbest)
#                2.5 %     97.5 %
#.sig01       0.000000  1.7042700
#.sigma       6.138221  7.5316389
#(Intercept) -2.230931 -0.2626198

# now compare adult female gonado-somatic index (GSI) across treatments and generations in warmed populations only.
# this analysis compares fisheries selection treatments to fisheries size selection controls (uniform size selection) in warmed populations only

# load and prepare data for warmed 'F' generation analysis
Dev_GSI_dataFH <- read.csv(file = "InputsZf/GSI_Dev_data_Warmed.csv")

Dev_GSI_dataFH$Treat <- as.factor(paste(Dev_GSI_dataFH$Temperature,Dev_GSI_dataFH$Selection, sep = "_"))

Dev_GSI_dataFH$Treat <- relevel(Dev_GSI_dataFH$Treat , ref="H_C")

# now for the model
GSIFMH <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_GSI_dataFH, REML = T)# 

summary(GSIFMH)

plot(GSIFMH) # check resid spread 

qqnorm(residuals(GSIFMH))

# now run the model selection for Deviance. 
# specify the full model
GSIFMH1a <- lmer(Deviance ~ Treat*Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_GSI_dataFH, REML = F) 
# specify the additive model
GSIFMH1b <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_GSI_dataFH, REML = F) 

# test the models using anova
anova(GSIFMH1a, GSIFMH1b)
#Data: Dev_GSI_dataFH1
#Models:
#GSIFMH1b: Deviance ~ Treat + Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#GSIFMH1a: Deviance ~ Treat * Generation + (1 | Population) + (1 | Tank) + (1 | Generation)
#         npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
#GSIFMH1b    8 2320.1 2351.3 -1152.0   2304.1                     
#GSIFMH1a   10 2323.0 2362.0 -1151.5   2303.0 1.1012  2     0.5766

# now fit the supported model
GSIFMH1best <- lmer(Deviance ~ Treat + Generation  + (1|Population) + (1|Tank) + (1|Generation), data = Dev_GSI_dataFH, REML = T)# 

plot(allEffects(GSIFMH1best)) # gives good visual of parameter estimates

summary(GSIFMH1best)
#Fixed effects:
#                 Estimate Std. Error       df t value Pr(>|t|)  
#(Intercept)        0.9060     0.8812 361.0000   1.028   0.3045  
#TreatH_G          -1.4353     0.7090 361.0000  -2.024   0.0437 *
#TreatH_S          -0.4569     0.7563 361.0000  -0.604   0.5461  
#Generation        -0.3428     0.2722 361.0000  -1.260   0.2086  
 

confint(GSIFMH1best)
#                      2.5 %      97.5 %
#.sig01       0.0000000  0.88362368
#.sig02       0.0000000  0.78430559
#.sig03       0.0000000  1.11256470
#.sigma       5.2936600  6.12070716
#(Intercept) -0.8126065  2.62353115
#TreatH_G    -2.8225181 -0.05230497
#TreatH_S    -1.9324885  1.01849870
#Generation  -0.8751115  0.18913038

#### Yield per recruit analysis ####
# functions
# modified fishmethods function
ypr_fishmethods <- function (age = NULL, wgt = NULL, partial = NULL, M = NULL, 
                             plus = FALSE, oldest = NULL, maxF = 2, incrF = 0.001, graph = TRUE) 
{
  if (is.null(age)) 
    stop("age vectoris missing")
  if (is.null(age)) 
    stop("wgt vector is missing")
  if (is.null(partial)) 
    stop("partial recruitment vector is missing")
  if (is.null(M)) 
    stop("M value or vector is missing")
  if (plus == TRUE & is.null(oldest)) 
    stop("oldest must be specified for plus group calculation.")
  if (any(length(age) != c(length(age), length(wgt), length(partial)))) 
    stop("Length of vectors unequal")
  if (length(M) == 1) 
    M <- rep(M, length(age))
  data <- as.data.frame(cbind(age, wgt, partial, M))
  YPR <- as.data.frame(cbind(rep(NA, ceiling(maxF/incrF) + 
                                   1), rep(NA, ceiling(maxF/incrF) + 1)))
  names(YPR) <- c("F", "YPR")
  if (plus == TRUE) {
    len <- oldest - min(data$age) + 1
    if (oldest > max(data$age)) {
      pdata <- data[rep(length(data$age), times = oldest - 
                          data$age[length(data$age)]), ]
      pdata$age <- seq(max(data$age) + 1, oldest, 1)
      data <- rbind(data, pdata)
    }
  }
  if (plus == FALSE) 
    len <- max(data$age) - min(data$age) + 1
  F <- 0
  for (i in 1:length(YPR$F)) {
    data$S <- exp(-data$partial * F - data$M)
    data$psb[1] <- 1
    for (y in 2:length(data$psb)) {
      data$psb[y] <- data$S[y - 1]
    }
    data$psb <- cumprod(data$psb)
    data$YPR <- ((data$partial * F)/(data$partial * F + 
                                       data$M)) * (1 - exp(-data$partial * F - data$M)) * 
      data$psb * data$wgt
    YPR$YPR[i] <- sum(data$YPR)
    YPR$F[i] <- F
    F <- F + incrF
  }
  Ymax <- max(YPR$YPR)
  Fmax <- YPR$F[which(YPR$YPR == max(YPR$YPR))]
  s10 <- ((YPR$YPR[2] - YPR$YPR[1])/(YPR$F[2] - YPR$F[1])) * 
    0.1
  F10 <- Fmax/2
  df <- F10/2
  ok <- 0
  fuzz <- 1e-04
  while (ok == 0) {
    data$s1 <- exp(-data$partial * F10 - data$M)
    data$s1p <- 1
    for (y in 2:length(data$s1p)) {
      data$s1p[y] <- data$s1[y - 1]
    }
    data$s1p <- cumprod(data$s1p)
    d1 <- sum(((data$partial * F10)/(data$partial * F10 + 
                                       data$M)) * (1 - exp(-data$partial * F10 - data$M)) * 
                data$s1p * data$wgt)
    data$s2 <- exp(-data$partial * (F10 + 1e-04) - data$M)
    data$s2p <- 1
    for (y in 2:length(data$s2p)) {
      data$s2p[y] <- data$s2[y - 1]
    }
    data$s2p <- cumprod(data$s2p)
    d2 <- sum(((data$partial * (F10 + 1e-04))/(data$partial * 
                                                 (F10 + 1e-04) + data$M)) * (1 - exp(-data$partial * 
                                                                                       (F10 + 1e-04) - data$M)) * data$s2p * data$wgt)
    slope <- (d2 - d1)/((F10 + 1e-04) - F10)
    if (abs(s10 - slope) <= fuzz) 
      ok <- 1
    if (ok == 0) {
      if (slope > s10) 
        F10 <- F10 + df
      if (slope < s10) 
        F10 <- F10 - df
      df <- df/2
    }
    Y10 <- d1
  }
  ans <- NULL
  ans <- matrix(NA, 2L, 2L)
  ans <- rbind(cbind(F10, Y10), cbind(Fmax, Ymax))
  dimnames(ans) <- list(cbind("F0.10", "Fmax"), c("F", "Yield_Per_Recruit"))
  outpt <- list(ans, YPR)
  names(outpt) <- c("Reference_Points", "F_vs_YPR")
  if (graph == TRUE) 
    plot(YPR[, 2] ~ YPR[, 1], ylab = "Yield-Per-Recruit", 
         xlab = "F", type = "l")
  return(outpt)
}

# Yield per recruit calculation function
ypr_func <- function(x, selCurve){
  x <- x
  ages <- x$Week
  ages2 <- x$Age
  weights <- x$weight
  selCurve = logist(inL50 = 10, delta = 4, depend = ages)
  plot(ages, selCurve, lwd = 2, type = "l", xlab = "age", ylab = "selectivity")
  yp <- ypr_fishmethods (age = ages, wgt = weights, partial = selCurve, M = 0.1, maxF = 5, incrF = 0.1)
  output <- data.frame(max_yield =  yp$Reference_Points[2,2], max_yield_f = yp$Reference_Points[2,1], F_point1 = yp$Reference_Points[1,2], F_point1_F = yp$Reference_Points[1,1], yield_F4 = yp$F_vs_YPR[41,2], yield_F2 = yp$F_vs_YPR[21,2])
  return(output)
}

# load data
# we want the F5 generation for temperature treatments, and the CG2 generation for fisheries treatments.
# CG2 is required as this is the first generation that we measured weight in that was not impacted by our size selection treatments.

weight_data_fishing_CG2 <- read.csv(file = "InputsZf/Weights_fishing_CG2.csv")
weight_data_warming_F5 <- read.csv(file = "InputsZf/Weights_warming_F5.csv")

Weight_dat <- bind_rows(weight_data_warming_F5, weight_data_fishing_CG2)

# estimate mean weight per week for each population (and each generation)

Weight_summary <- Weight_dat %>% dplyr::group_by(Population, Tank, Selection, Temperature, Generation, Age, Week) %>% dplyr::summarise(weight = mean(Weight)) %>% drop_na()

# now predict yield 

#select age at 50% selectivity - change this value later for the exercise
#Normally we would know length and not age of 50% selectivity, but here for simplicity we assume that we calcuated age from length already 
ag50 = 10
#delta is the age difference between 50 and 95% (or 5% and 50%) selectivity
delta = 2
#now we will use 'logist' function from MQMF package to calculate selectivity at each age. First check how the function works 
#?logist

# create a selectivity curve
selCurve = logist(inL50 = ag50, delta = delta, depend = c(min(Weight_summary$Week):max(Weight_summary$Week)))
#and plot the curve
plot(c(min(Weight_summary$Week):max(Weight_summary$Week)), selCurve, lwd = 2, type = "l", xlab = "age", ylab = "selectivity")

# now calculate YPR for each population in each generation
ypr_pops <- Weight_summary %>% group_by(Population, Generation, Temperature, Selection) %>% do(ypr_func(., selCurve))
# summarise across treatments
ypr_pops_summary <- ypr_pops %>% group_by(Generation, Temperature, Selection) %>% summarise(mean_max = mean(max_yield), sd_max = sd(max_yield) ,mean_point1 = mean(F_point1), sd_point1 = sd(F_point1) , mean_F4 = mean(yield_F4), sd_F4 = sd(yield_F4), mean_F2 = mean(yield_F2), sd_F2 = sd(yield_F2))
ypr_pops_summary

# A tibble: 5 × 11
# Groups:   Generation, Temperature [3]
#  Generation Temperature Selection mean_max  sd_max mean_point1 sd_point1 mean_F4  sd_F4 mean_F2   sd_F2
#.      <int> <chr>       <chr>        <dbl>   <dbl>       <dbl>     <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
#1          5 H           C            0.417 0.0804        0.375    0.0770   0.414 0.0771   0.409 0.0829 
#2          5 L           C            0.487 0.0377        0.446    0.0400   0.477 0.0272   0.483 0.0404 
#3          7 L           C            0.430 0.0768        0.396    0.0725   0.408 0.0612   0.426 0.0725 
#4          7 L           G            0.516 0.00675       0.471    0.0211   0.493 0.0143   0.507 0.00302
#5          7 L           S            0.384 0.0483        0.347    0.0338   0.376 0.0532   0.382 0.0475 


