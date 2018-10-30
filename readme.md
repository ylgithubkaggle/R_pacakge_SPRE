SPRE
================

Download and install the package from github
--------------------------------------------

The whole process can be done in R. First check the availability of "devtools":

    if(!require("devtools"))install.packages("devtools")
    library(devtools)

After installing devtools, the SPRE package can be installed simply by running:

    install_github("ylgithubkaggle/R_pacakge_SPRE",subdir="SPRE")

Finally, everytime we need to use the package, we can load the package by running:

    library(SPRE)

View tutorial and help files
----------------------------

A vignette will tell you how to do the analysis. To view the vignette, we can do:

    vignette("SPRE")

A help file can be loaded by running:

    help(package="SPRE")

We can load the original data by running:

    data("HEAT_stat")

SPRE is the function that does the major analysis, to get help with the function, we can do

    SPRE

We can perform SPRE analyses by running:

    modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
    plot_residuals(modelfit)
    stability(modelfit)
    plot_weibull(modelfit)

The outputs can also be found in the vignette.
