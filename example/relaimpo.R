library(relaimpo)
#This R script conducts relative importance decomposition to assign shares of importance of each variable and plot the results.

data=read.csv('./example/relative_contribution_ILS_Err_Intro.csv')
bootimpo.result <- boot.relimp(data, b = 100,
                    type = c("lmg", "last", "first", "pratt"),
                    rank = TRUE, diff = TRUE, rela = TRUE)

plot(bootimpo.result)