library(testthat)
library(MANOVA.RM)
context("multRM")

test_that("multRM equals RM",{
  data(EEG)
  t1 <- multRM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", 
               within = c("region", "feature"), iter = 1, CPU =1)
  t2 <- RM(resp ~ diagnosis * region * feature, data = EEG, subject = "id", 
               no.subf = 2, iter = 1, CPU = 1)
  expect_equal(t1$WTS, t2$WTS)
})


test_that("missing values", {
  data(EEG)
  EEG2 <- EEG[-5, ]
  expect_error(multRM(resp ~ sex * feature * region, data = EEG2,
                      within = c("feature", "region"),
                  subject = "id", iter = 1, CPU = 1))
})
