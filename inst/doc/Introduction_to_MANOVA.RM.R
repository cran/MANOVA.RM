## ------------------------------------------------------------------------
library(MANOVA.RM)
data(o2cons)

## ------------------------------------------------------------------------
head(o2cons)

## ------------------------------------------------------------------------
set.seed(1234)
model1 <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons, 
             subject = "Subject", no.subf = 2, iter = 10000, resampling = "Perm", CPU = 1)
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
set.seed(987)
EEG_MANOVA <- MANOVA(resp ~ sex * diagnosis, 
                     data = EEG, subject = "id", resampling = "paramBS", 
                     iter = 1000,  alpha = 0.01, CPU = 1)
summary(EEG_MANOVA)

## ------------------------------------------------------------------------
GUI.MANOVA()

