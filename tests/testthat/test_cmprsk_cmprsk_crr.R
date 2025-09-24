skip_if_not_installed("cmprsk")
skip_if_not_installed("mlr3proba")

test_that("autotest", {
  message("\n =====> autotest")

  learner = lrn("cmprsk.crr", maxiter = 99, gtol = 1e-3)
  print(learner)
  #expect_learner(learner)
  result =run_autotest(  
    learner,
    N = 42,
    check_replicable = FALSE
  )
  expect_true(result, info = result$error)
  message(" autotest ENDS")
})

