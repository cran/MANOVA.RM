## ------------------------------------------------------------------------
library(MANOVA.RM)
data(o2cons)

## ------------------------------------------------------------------------
head(o2cons)

## ------------------------------------------------------------------------
model1 <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons, 
             subject = "Subject", no.subf = 2, iter = 10000, 
             resampling = "Perm", CPU = 1, seed = 1234)
summary(model1)

## ------------------------------------------------------------------------
data(EEG)
set.seed(987)
EEG_model <- RM(resp ~ sex * diagnosis * feature * region, 
                     data = EEG, subject = "id", no.subf = 2, resampling = "WildBS",
                     iter = 1000,  alpha = 0.01, CPU = 1)
summary(EEG_model)

## ------------------------------------------------------------------------
plot(EEG_model, factor = "sex", main = "Effect of sex on EEG values")
plot(EEG_model, factor = "sex:diagnosis", legendpos = "topleft", col = c(4, 2))
plot(EEG_model, factor = "sex:diagnosis:feature", legendpos = "center")

## ------------------------------------------------------------------------
data(EEG)
EEG_MANOVA <- MANOVA(resp ~ sex * diagnosis, 
                     data = EEG, subject = "id", resampling = "paramBS", 
                     iter = 1000,  alpha = 0.01, CPU = 1, seed = 987)
summary(EEG_MANOVA)

## ------------------------------------------------------------------------
library(HSAUR)
data(water)
water$subject <- 1:dim(water)[1]
library(tidyr)
data_long <- gather(water, measurement, response, mortality:hardness, factor_key=TRUE)
head(data_long)
test <- MANOVA(response ~ location, data = data_long, subject = "subject", iter = 1000, resampling = "paramBS", CPU = 1, seed = 123)
summary(test)
cr <- conf.reg(test)
cr

## ------------------------------------------------------------------------
plot(cr, col = 2, lty = 2, xlab = "Difference in mortality", ylab ="Difference in water hardness")

## ------------------------------------------------------------------------
GUI.MANOVA()

