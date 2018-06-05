
context("test_fit_diagnostics")

test_that("analyze fit diagnostics", {
    pdf("test1.pdf")
    retina_A <- main_3hgbqg("3hgbqg" %>% get_path_to_test_retina_folder)
    retina_B <- main_40oik5("40oik5" %>% get_path_to_test_retina_folder)
    retina_C <- main_gv3igs("gv3igs" %>% get_path_to_test_retina_folder)
    
    retinaplot(retina_A)
    fit_plots(retina_A$fit_data1)
    fit_plots(retina_B$fit_data1)
    fit_plots(retina_C$fit_data1)
    dev.off()
    
})
