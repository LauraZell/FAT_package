## Created March 13, 2025
## Modified March
## Last modified March


# Developing R package Forecasted Average Treatment Method
# ----------------------------------
## SECTION 1: LOAD PACKAGES
# ----------------------------------
# Include all dependencies. Select "NO" when asked to install from source
install.packages(c("devtools","tidyverse","usethis","testthat","knitr","roxygen2"))

library(pacman)
pacman::p_load(devtools,tidyverse,usethis,testthat,knitr,roxygen2)

sessionInfo() # Check version information about R and OS
# I am using: R (version: 4.4.3), macOS Squioa 15.3.1, tidyverse (version: 2.0.0), devtools (version 2.4.5)

# ----------------------------------
## SECTION 2: CREATE PACKAGE
# ----------------------------------
#create_package("/Users/laurazell/Dropbox/04_Doktorat/research_stay/UCL/R_package/fat/fatpackage")

#document() #command + shift + d -> create and update documentation
#check() #command + shift + e -> build and check package


# ----------------------------------
## SECTION 3: CREATE FUNCTION
# ----------------------------------

# ----------------------------------
## SECTION 4: TESTING THE FUNCTION
# ----------------------------------
use_testthat()

