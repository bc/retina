context("count_to_rho")

test_that("count_to_rho converts raw counts per square micron to counts per square millimeter", {
expect_equal(count_to_rho(15, 25, 25), 24000)
expect_equal(count_to_rho(10, 10, 10), 1e+05)
expect_equal(count_to_rho(0, 10, 10),   0)
})




