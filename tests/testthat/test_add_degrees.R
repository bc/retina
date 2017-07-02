context("add_degrees")

test_that("add_degrees spins properly", {
expect_true(abs(add_degrees(0, 0) -  0)		 <= .00000001)
expect_true(abs(add_degrees(-90, 0) -  -90)	 <= .00000001)
expect_true(abs(add_degrees(-180, 0) -  -180)<= .00000001)
expect_true(abs(add_degrees(-180, 90) -  -90)<= .00000001)
expect_true(abs(add_degrees(180, 360) -  180)<= .00000001)
expect_true(abs(add_degrees(150, 40) -  -170)<= .00000001)
expect_true(abs(add_degrees(90, -100) -  -10)<= .00000001)


})
