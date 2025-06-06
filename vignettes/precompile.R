# Execute the code from the vignette
knitr::knit("vignettes/guide.Rmd.orig", output = "vignettes/guide.Rmd")
file.rename("figure/target-effects-plot-1.png", "vignettes/figure/target-effects-plot-1.png")
