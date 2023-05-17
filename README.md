# GMNL Package User Guide

This guide will provide you with the necessary steps to get started with the GMNL package. It will also guide you through data formatting, McFadden $R^2$, $\rho^2$ ASC, Robust/Standard errors, confidence, and prediction intervals in WTP space.

## Table of Contents

- [Getting Started](#getting-started)
- [Data Formatting](#data-formatting)
- [McFadden R^2](#mcfadden-r2)
- [Rho^2 ASC](#rho2-asc)
- [Robust/Standard Errors](#robuststandard-errors)
- [Confidence and Prediction Intervals in WTP Space](#confidence-and-prediction-intervals-in-wtp-space)

## Getting Started

To get started with the GMNL package, you'll need to download the package zip file. Make sure once you download it, you delete the README.md file from the folder you launch the R Package from to prevent any errors. It is assumed that you already have R and RStudio installed on your machine to use this package.

Once you've downloaded the zip file, you can open the GMNL package by following these steps:

1. Open the zip folder you downloaded from this page
2. Double click the `gmnl.Rproj` file 
3. This should open RStudio. Next, navigate to the Build tab in the top right, then click Install. 
4. The Console should indicate that your R session has restarted and automatically run `library(gmnl)` 
5. You are ready to use gmnl!

Your GMNL package is now ready for use!

## Data Formatting

The GMNL package requires your data to be in a specific format. If you have previously used the `mlogit.data()`, this will be the same but instead you will call `gmnl.data()`. `gmnl.data()` uses the same inputs, so if you have existing code you simply can replace the function name and it will work. 

```markdown
data("TravelMode", package = "AER") # import the dataset 

TM <- gmnl.data(TravelMode,choice = "choice", shape = "long", alt.levels = c("air", "train", "bus", "car"))

# For a wide dataset

data("Electricity", package = "mlogit") # import the dataset 2

Electr <- gmnl.data(Electricity, id.var = "id", choice = "choice", varying = 3:26, shape = "wide", sep = "")

``` 

## McFadden R^2

The McFadden R^2 is a measure of fit for the model. It is computed as follows:

```markdown
# Format out data to the gmnl data type 
NYC <- gmnl.data(NYC_short, 
                id.var = "ResponseId", 
                choice = "choice", 
                varying = 3:20, 
                shape = "wide" )

# Create a MIXL model                
mod1 <- gmnl(choice ~ price + crowd + wait + emissions , 
            data = NYC, 
            model = "mixl", 
            R=30, 
            panel = TRUE,
            ranp = c(price = "n", crowd = "n", wait = "n", emissions = "n")
            )
# Print out our model summary
summary(mod1)

```

