context("add_degrees")

test_that("add_degrees spins properly", {
    expect_true(abs(add_degrees(0, 0) - 0) <= 1e-08)
    expect_true(abs(add_degrees(-90, 0) - -90) <= 1e-08)
    expect_true(abs(add_degrees(-180, 0) - -180) <= 1e-08)
    expect_true(abs(add_degrees(-180, 90) - -90) <= 1e-08)
    expect_true(abs(add_degrees(180, 360) - 180) <= 1e-08)
    expect_true(abs(add_degrees(150, 40) - -170) <= 1e-08)
    expect_true(abs(add_degrees(90, -100) - -10) <= 1e-08)
    
    
})
