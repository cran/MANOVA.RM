library(testthat)
library(MANOVA.RM)
context("RM output")

test_that("example 1: 1 whole, 2 sub",{
  oxy <- RM(O2 ~ Group * Staphylococci * Time, data = o2cons, 
            subject = "Subject", no.subf = 2, iter = 100, CPU = 1)
  expect_equal(oxy$WTS[1], 11.167)
})


test_that("example 2: 2 whole, 2 sub", {
  EEG_model <- RM(resp ~ sex * diagnosis * feature * region,
                  data = EEG, subject = "id", no.subf = 2, resampling = "WildBS",
                  iter = 100, CPU = 1)
  expect_equal(EEG_model$WTS[1], 9.973)
})


test_that("missing interaction",{
  data(EEG)
  expect_error(RM(resp ~ sex + feature + region, data = EEG,
                 subject = "id", iter = 10, CPU = 1))
})

test_that("wrong interaction",{
  data(EEG)
  expect_error(RM(resp ~ sex + feature * region, data = EEG,
                  subject = "id", iter = 10, CPU = 1))
})

test_that("missing values", {
  data(EEG)
  EEG2 <- EEG[-5, ]
  expect_error(RM(resp ~ sex * feature * region, data = EEG2, no.subf = 2,
                  subject = "id", iter = 1, CPU = 1))
})
