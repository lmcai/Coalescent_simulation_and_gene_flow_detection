library(relaimpo)
#This R script conducts relative importance decomposition to assign shares of importance of each variable and plot the results.

x=read.csv('./example/relative_contribution_ILS_Err_Intro.csv')
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
#calculate relative importance
calc.relimp(data, type = "lmg",rela=TRUE)

#bootstrap for confidence intervals
bootimpo.result <- boot.relimp(Y_BSsum~X1_ILS+X2_Err+X3_Hyb,data=data, b = 100,
                    type = c("lmg", "last", "first", "pratt"),
                    rank = TRUE, diff = TRUE, rela = TRUE)

booteval.relimp(bootimpo.result,lev=0.9,nodiff=TRUE)

## Plot
plot(booteval.relimp(bootimpo.result))

##Adding interaction terms in the regression model 

linmod <- lm(Y_BSsum ~ log(X1_ILS)+log(X2_Err)+log(X3_Hyb)+log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * X3_Hyb)+log(X1_ILS * X3_Hyb*X2_Err), data = data)
summary(linmod)

#Call:
#lm(formula = Y_BSsum ~ log(X1_ILS) + log(X2_Err) + log(X3_Hyb) + 
#    log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * 
#    X3_Hyb) + log(X1_ILS * X3_Hyb * X2_Err), data = data)

#Residuals:
#    Min      1Q  Median      3Q     Max 
#-466.09  -84.64   14.36  104.33  600.90 

#Coefficients: (4 not defined because of singularities)
#                              Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                     482.05     117.30   4.109 0.000124 ***
#log(X1_ILS)                     -59.49      33.67  -1.767 0.082412 .  
#log(X2_Err)                     106.12      12.64   8.394 1.18e-11 ***
#log(X3_Hyb)                      66.92      24.25   2.759 0.007702 ** 
#log(X1_ILS * X2_Err)                NA         NA      NA       NA    
#log(X1_ILS * X3_Hyb)                NA         NA      NA       NA    
#log(X2_Err * X3_Hyb)                NA         NA      NA       NA    
#log(X1_ILS * X3_Hyb * X2_Err)       NA         NA      NA       NA    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 188 on 59 degrees of freedom
#Multiple R-squared:  0.6694,	Adjusted R-squared:  0.6526 
#F-statistic: 39.83 on 3 and 59 DF,  p-value: 3.363e-14
