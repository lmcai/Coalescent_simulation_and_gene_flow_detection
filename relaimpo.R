library(relaimpo)
#This R script conducts relative importance decomposition to assign shares of importance of each variable and plot the results.

x=read.csv('./example/relative_contribution_ILS_Err_Intro.csv')
data=x[,c('Y_BSsum','X1_ILS','X2_Err','X3_Hyb')]
bootimpo.result <- boot.relimp(data, b = 100,
                    type = c("lmg", "last", "first", "pratt"),
                    rank = TRUE, diff = TRUE, rela = TRUE)

## Plot
plot(booteval.relimp(bootimpo.result))
