# GMNL Package User Guide

This guide will provide you with the necessary steps to get started with the GMNL package. It will also guide you through data formatting, McFadden $R^2$, $\rho^2$ ASC, Robust/Standard errors, confidence, and prediction intervals in WTP space.

## Table of Contents

- [Getting Started](#getting-started)
- [Data Formatting](#data-formatting)
- [McFadden $R^2$](#mcfadden-r2)
- [$\rho^2$ ASC](#rho2-asc)
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

## McFadden $R^2$

The McFadden $R^2$ is a measure of fit for the model. It is computed as follows:

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

It is possible to simply return the McFadden values by calling the respective McFadden $R^2$ and McFadden $R^_{Adj}$ functions on the returned model object.

```markdown
mcfadden.gmnl(mod1)
# 0.216984

mcfadden_adj.gmnl(mod1)
# 0.2145559 
```

We can also return the value inside a new GMNL object and then access it.

```
mod1_mcf <- mcfadden.gmnl(mod1, returnObj = TRUE)
mod1_mcf$mcfadden
#[1] 0.216984

mod1_mcfA <- mcfadden_adj.gmnl(mod1, returnObj = TRUE)
mod1_mcfA$mcfaddenadj
#[1] 0.2145559
```


## $\rho^2$ ASC

Similar to McFadden, $\rho^2$ allows us to see the performance of our model in terms of capturing unobserved preference traits of our population. We can run a model on the same data as the above where we include alternative specific constants in our model to see if we are able to better predict preference. We have also included the feature of simply returning the $\rho^2$ or ${\rho^2}_{Adj}$ value by directly calling the functions associated with calculating the respective value.

$\rho^{2}_{ASC} = 1 - \frac{L(\hat{\beta})}{L(ASC)}$

${\rho^{2}_{ASC}}_{adj} = 1 - \frac{L(\hat{\beta}) - K}{L(ASC)}$

This statistic functions in a similar manner to McFadden’s where the larger the value, the better the model for the particular set of data. The difference between the two statistics is that $L(ASC)$ is not a null model. We use $L(ASC)$ to represent the model where the choice probabilities are constant across all alternatives and no predictors are considered. In terms of the name of $L(ASC)$ being used instead of something more intuitive like $L(ll_1)$, this model is useful to help compare a model with alternative specific constants to one without. Although this statistic is displayed for models both with and without alternative specific constants, it helps evaluate if the use of alternative specific constants allow for the model to capture preference which is not explained by an explanatory variable included in the model. Accounting for unobserved factors can improve the model fit. For the adjusted version of this model, it is similar to McFadden's where it penalizes for the amount of explanatory variables inclued in the model. 


```
mod2 <- gmnl(choice ~ price + crowd + wait  | 1 | 0 , 
            data = NYC, 
            model = "mixl", 
            R=30, 
            panel = TRUE,
            ranp = c(price = "n", crowd = "n", wait = "n"))

summary(mod2)

# Return the Rho^2 and Rho^2 adjusted values

rho2.gmnl(mod2)
# 0.2115554

rho2_adj.gmnl(mod2)
# 0.209064

# Return and store GMNL object with Rho Squared

mod2_rho <- rho2.gmnl(mod2, returnObj = TRUE)
mod2_rho$rho2
#[1] 0.2115554

mod2_rhoA <- rho2_adj.gmnl(mod2, returnObj = TRUE)
mod2_rhoA$rho2adj
#[1] 0.209064
```
## Robust/Standard Errors

The robust standard errors GMNL includes a covariance matrix of the standard errors
for the computed model. The matrix computed is symmetrical over the diagonal and is
$n$ x $n$ in size where $n$ is the number of coefficients in the model. The values on the
diagonal are the variances of the respective random variables, while the variables off
the diagonal are the covariance of each pair of variables.

```
# We continue with the model we had previously calculated 
compute_clustered_standard_error.gmnl(mod1)
```

## Confidence and Prediction Intervals in WTP Space

The $Delta_{SMG}$ method proposed by [Scaccia, Marcucci, and Gatta (2023)](https://www.sciencedirect.com/science/article/pii/S0191261522001886) address the
time intensive method of computing the KR procedure as it involves Monte Carlo
simulations in many dimensions. Their approach builds on KR and another method
proposed by Bliemer and Rose (2013), $Delta_{BR}$. $Delta_{BR}$ avoids much of the
simulation required for KR through deriving partial analytical expressions for the
standard errors. $Delta_{SMG}$ combines aspects from both KR and $Delta_{BR}$ such that it
produces the whole WTP simulated sampling distribution and is partly analytical
which reduces the computational burden.

Let us consider that WTP is defined by $w_k$ and is no longer a random variable, but
instead defined by the following ratio:

$w_k(z) = - \frac{\beta_k(\theta_k | z_k)}{\beta_c(\theta_c|z_c)}$

Such that theta is a fixed parameter and z is a random variable that accounts for
heterogeneity. We also enforce that $\beta_k$ and $\beta_c$ are expressed in parameter free
form, meaning they follow parameter free distributions such as the standard normal or
standard uniform distribution.

Next, calculate parameter vector estimates and estimated variance covariance matrix of
thetas. A large number of pseudo random, such as Halton draws, are obtained for z
from the previously selected standard distribution where ${z_1}^r$ where $r=1 \dots R$ denote the draws. For each $r$ from the set of draws, compute $\beta_k(\theta_k,z_k)$ and then compute the willingness to pay estimate at $z_r$ where:

$\hat{w}_k(z^{(r)}) = - \frac{\beta_k(\hat{\theta}_k , z_k^{(r)})}{\beta_c(\hat{\theta}_c,z_c^{(r)})}$

The resulting vector of $w_k$ can be used to get different values of the confidence interval such as average WTP, standard error,confidence intervals, prediction standard error, and estimated prediction intervals.

For this section, methods for both Krinsky Robb and $Delta_{SMG}$ as functions to compute the confidence and prediction intervals of WTP space. 

For the Krinsky Robb method :
```
# Import the Electricity dataset from the mlogit package
data("Electricity", package = "mlogit")

Electr <- gmnl.data(Electricity, 
                    id.var = "id", 
                    choice = "choice",
                    varying = 3:26, 
                    shape = "wide", 
                    sep = "")

# Estimate a MIXL model with random parameters
Elec.gmnl <- gmnl(choice ~ pf + cl + loc + wk + tod + seas| 0, 
                  data = Electr, 
                  subset = 1:3000, 
                  model = 'mixl', 
                  R = 10, 
                  panel = TRUE, 
                  ranp = c(cl = "n", loc = "n", wk = "n", 
                           tod = "n", seas = "n"), 
                  correlation = FALSE)

# Estimate WTP and CIs using Krinsky and Robb method
KRint=KRCI(Elec.gmnl,
           Q = 5000,
           B = 2000,
           kcoe = c(2:6),
           costcoe = 1,
           level = 0.95,
           qw = 1.5,
           tpw = 0.1)

# Print out summary of the estimates
summary.KRCI(KRint)
```

For the new proposed $Delta_SMG$ method:

```
#Compute the confidence intervals
Dint=DeltaCI(Elec.gmnl,
             Q = 5000,
             kcoe = c(2:6),
             costcoe = 1,
             level = 0.95,
             qw = 1.5,
             tpw = 0.1)

# Call the summary function for Delta method
summary.Delta(Dint)
```

