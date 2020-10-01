library(relaimpo)
#This R script conducts relative importance decomposition to assign shares of importance of each variable and plot the results.

x=read.csv('./example/relative_contribution_ILS_Err_Intro.csv')
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
bootimpo.result <- boot.relimp(data, b = 100,
                    type = c("lmg", "last", "first", "pratt"),
                    rank = TRUE, diff = TRUE, rela = TRUE)

## Plot
plot(booteval.relimp(bootimpo.result))

##Adding interaction terms in the regression model 
linmod <- lm(Y_BSsum ~ X1_ILS+X2_Err+X3_Hyb+X1_ILS*X2_Err+X1_ILS*X3_Hyb+X2_Err*X3_Hyb, data = data)

summary(linmod)

#Call:
#lm(formula = Y_BSsum ~ X1_ILS + X2_Err + X3_Hyb + X1_ILS * X2_Err + 
#    X1_ILS * X3_Hyb + X2_Err * X3_Hyb, data = data)

#Residuals:
#    Min      1Q  Median      3Q     Max 
#-393.99  -83.40   -6.18   94.61  756.73 

#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)   
#(Intercept)    1.945e+02  7.873e+01   2.471  0.01656 * 
#X1_ILS         2.775e+00  1.586e+00   1.750  0.08562 . 
#X2_Err         9.297e+00  2.820e+00   3.297  0.00170 **
#X3_Hyb         8.737e+00  3.074e+00   2.842  0.00623 **
#X1_ILS:X2_Err -2.869e-02  3.662e-02  -0.783  0.43676   
#X1_ILS:X3_Hyb -8.002e-04  5.376e-02  -0.015  0.98818   
#X2_Err:X3_Hyb -1.023e-01  4.156e-02  -2.462  0.01693 * 
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 202.1 on 56 degrees of freedom
#Multiple R-squared:  0.6374,	Adjusted R-squared:  0.5985 
#F-statistic: 16.41 on 6 and 56 DF,  p-value: 8.486e-11

linmod <- lm(Y_BSsum ~ log(X1_ILS)+log(X2_Err)+log(X3_Hyb)+log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * X3_Hyb)+log(X1_ILS * X3_Hyb*X2_Err), data = data)
summary(linmod)

Call:
lm(formula = Y_BSsum ~ log(X1_ILS) + log(X2_Err) + log(X3_Hyb) + 
    log(X1_ILS * X2_Err) + log(X1_ILS * X3_Hyb) + log(X2_Err * 
    X3_Hyb) + log(X1_ILS * X3_Hyb * X2_Err), data = data)

Residuals:
    Min      1Q  Median      3Q     Max 
-466.09  -84.64   14.36  104.33  600.90 

Coefficients: (4 not defined because of singularities)
                              Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     482.05     117.30   4.109 0.000124 ***
log(X1_ILS)                     -59.49      33.67  -1.767 0.082412 .  
log(X2_Err)                     106.12      12.64   8.394 1.18e-11 ***
log(X3_Hyb)                      66.92      24.25   2.759 0.007702 ** 
log(X1_ILS * X2_Err)                NA         NA      NA       NA    
log(X1_ILS * X3_Hyb)                NA         NA      NA       NA    
log(X2_Err * X3_Hyb)                NA         NA      NA       NA    
log(X1_ILS * X3_Hyb * X2_Err)       NA         NA      NA       NA    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 188 on 59 degrees of freedom
Multiple R-squared:  0.6694,	Adjusted R-squared:  0.6526 
F-statistic: 39.83 on 3 and 59 DF,  p-value: 3.363e-14