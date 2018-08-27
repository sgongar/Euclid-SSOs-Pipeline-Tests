# Read CSV into R
SSOsDataFrame <- read.csv(file="ssos_2.csv", header=TRUE, sep=",")

library(MASS)
library(ggplot2)
library(boot)

# cv_MSE_k10 <- rep(NA,10)
# 
# for (i in 1:10) {
#   modelo <- glm(B_IMAGE ~ poly(MAG_AUTO, i, raw=TRUE), data = SSOsDataFrame)
#   set.seed(17)
#   cv_MSE_k10[i] <- cv.glm(data = SSOsDataFrame, glmfit = modelo, K = 10)$delta[1]
# }
# p4 <- ggplot(data = data.frame(polinomio = 1:10, cv_MSE = cv_MSE_k10),
#              aes(x = polinomio, y = cv_MSE)) +
#       geom_point(colour = c("firebrick3")) +
#       geom_path()
# p4 <- p4 + theme(panel.grid.major = element_line(colour  =  'gray90'))
# p4 <- p4 + theme(plot.title = element_text(face  =  'bold'))
# p4 <- p4 + theme(panel.background = element_rect(fill  =  'gray98'))
# p4 <- p4 + labs(title  =  'Test Error ~ Grado del polinomio')
# p4 <- p4 + scale_x_continuous(breaks = 1:10)
# print(p4)


out <- ggplot(data = SSOsDataFrame, aes(x = MAG_AUTO, y = B_IMAGE)) +
              geom_point(color = "grey30", alpha = 0.3) + 
              geom_smooth(method = "lm", formula = y ~ poly(x, 5, raw=TRUE), color = "red") +
              labs(title = "Polinomio de grado 1: wage ~ age") +
              theme_bw() +
              theme(plot.title = element_text(hjust = 0.5))

# print(out)

modelo_poli5 <- lm(B_IMAGE ~ poly(MAG_AUTO, 5, raw=TRUE), data = SSOsDataFrame)
# print(summary(modelo_poli5))
# print(modelo_poli5$coefficients)

mpi <- cbind(SSOsDataFrame, predict(modelo_poli5, interval = "prediction"))
out_2 <- ggplot(mpi, aes(x = SSOsDataFrame$MAG_AUTO)) +
                geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "gray", alpha = 0.5) +
                geom_point(aes(y = SSOsDataFrame$B_IMAGE)) +
                geom_line(aes(y = fit), colour = "blue", size = 1)
# print(summary(mpi$fit))
# print(out_2)

out_3 <- predict(modelo_poli5, se.fit = TRUE, scale = NULL, df = Inf,
                 interval = c("prediction"), level = 0.95, type = c("response"))
# print(summary(out_3))
out_matrix <- out_3$fit
fit_list <- out_matrix[, "fit"]
lwr_list <- out_matrix[, "lwr"]
upr_list <- out_matrix[, "upr"]
print(max(fit_list))


