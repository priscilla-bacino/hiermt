#
#
#

# set.seed(1)
# n <- 50
#
# df <- data.frame(y1 = rnorm(n),
#                  y2 = rnorm(n),
#                  x = sample(rep(0:1,c(n/2,n/2))))
#
# hiermt(formula = cbind(y1,y2)~x,
#        data = df,
#        global_test = "bonferroni",
#        alpha = 0.05)
#
# hiermt(formula = .~x,
#        data = df,
#        global_test = "bonferroni",
#        alpha = 0.05)

test_that("this is a placeholder for actual tests", {
  skip("Not yet implemented")
})

# test_that("hiermt works as expected", {
#   expect_equal(class(hiermt(formula = .~x,
#                       data = df,
#                       global_test = "bonferroni",
#                       alpha = 0.05)), 'NULL')
# })

