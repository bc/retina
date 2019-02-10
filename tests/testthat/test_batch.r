
context("batch processing lizards")

lapply_testthat <- function(vec_of_strings, FUN) {
    lapply(vec_of_strings, function(e) {
        context(e)
        test_that(e, {
            FUN(e)
        })
    })
}
pdf("out.pdf")
q <- lapply_testthat(c("3hgbqg", "3hgbqg_direct_csv_coords", "8pm223", "8pm223_direct_csv_coords", 
    "bjii4t", "jhvbwt", "40oik5", "gv3igs", "raxp91"), run_diagram_retina_folder)
dev.off()