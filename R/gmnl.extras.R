utils::globalVariables(c("choice", "elast1", "elast2", "pred_prob", "prob_chg"))

# Prints out coefficients from model
coef.gmnl <- function(object, ...) {
  round(summary(object)$CoefTable, 3)
}

#####################################################################

#' Cumulative Density Function for Mixture of Normal Distributions
#'
#' @param x value at which CDF needs to be computed
#' @param w vector of weights for the normal densities in the mixture
#' @param u vector of means for the normal densities in the mixture
#' @param s vector of std dev for the normal densities in the mixture
#'
#' @return cumulative density function of x
#' @export
Fa <- function(x, w, u, s) {
  sum(w * pnorm(x, mean = u, sd = s))
}

# Compute quantiles for a mixture of normal distributions

#' Compute Quantiles for Mixture of Normal Distributions
#'
#' @param p CDF corresponding to the quantile to be computed
#' @param w vector of weights for the normal densities in the mixture
#' @param u vector of means for the normal densities in the mixture
#' @param s vector of std dev for the normal densities in the mixture
#' @param br vector of lower and upper bounds of the the interval in which the quantile is to be searched for
#'
#' @return the quantile corresponding to the cumulative density function p
#' @export
F_inv <- function(p, w, u, s, br = c(-1000, 1000)) {
  G = function(x)
    Fa(x, w, u, s) - p
  return(uniroot(G, br)$root)
}

#####################################################################


#' Computes Variety of Confidence Interval Data
#'
#' @param res an object of class glmn
#' @param Q number of random draws from standard distributions (normal or uniform)
#' @param kcoe vector of positions (in the res$coe vector) of the parameters for which wtp is to be computed
#' @param costcoe position (in the res$coe vector) of the cost parameter
#' @param level vector of confidence level to compute confidence and prediction intervals (i.e. 0.9, 0.95…)
#' @param qw value (the same for all parameter for which wtp is to be computed) at which upper tail probability
#'           for wtp distribution is to be computed
#' @param tpw probability (the same for all parameters for which wtp is to be computed) at which quantiles
#'            leaving that probability in the upper tail of wtp distribution are to be computed
#'
#' @return vector of : "WTP", "se.mean.WTP","pred.err.WTP","CI.mean.WTP","PI.WTP","quant.WTP","tail.prob.WTP"
#' WTP : Willingness to pay point estimate
#' se.mean.WTP : Willingness to pay point estimate
#' pred.err.WTP : Standard error for the mean WTP
#' CI.mean.WTP : Prediction error for WTP
#' PI.WTP : Prediction interval for WTP
#' quant.WTP : Quantile such that tpw% of the population has a larger WTP
#' tail.prob.WTP : Probability of observing a WTP larger than qw
#'
#' @export
DeltaCI <- function(res, Q, kcoe, costcoe, level, qw, tpw) {
  wtp_names <- res$wtp_names
  ncoe = length(kcoe) + length(costcoe)
  coe = c()
  rcoe = c()

  # Correcting for negative estimates of standard deviation
  # for random coefficients (sometimes happening with mlogit,
  # not sure this is required for gmnl)
  rp = length(res$ranp)
  for (i in (1:rp)) {
    res$coe[paste0("sd", sep = '.', names(res$ranp)[i])] = abs(res$coe[paste0("sd", sep =
                                                                                '.', names(res$ranp)[i])])
  }

  # Determins the type of distribution for paramaters of interest
  # and creates an ordered vector with coefficients (or means if random) of the attributes for which
  # wtp needs to be calculated, the cost coefficient (or mean if random)
  # and standard deviations of random coefficients
  typepar = rep(0, ncoe)
  for (i in (1:length(kcoe)))
  {
    a = names(res$coe)[kcoe[i]]
    if (sum(a == names(res$ranp)) == 0) {
      typepar[i] = "fixed"
      coe = c(coe, res$coe[a])
    } else if (res$ranp[a] == "ln") {
      typepar[i] = "ln"
      coe = c(coe, res$coe[a])
      rcoe = c(rcoe, res$coe[paste0("sd", sep = ".", a)])
    } else if (res$ranp[a] == "n") {
      typepar[i] = "n"
      coe = c(coe, res$coe[a])
      rcoe = c(rcoe, res$coe[paste0("sd", sep = ".", a)])
    }
  }
  a = names(res$coe)[costcoe]
  if (sum(a == names(res$ranp)) == 0) {
    typepar[length(kcoe) + 1] = "fixed"
    coe = c(coe, res$coe[a])
  } else if (res$ranp[a] == "ln") {
    typepar[length(kcoe) + 1] = "ln"
    coe = c(coe, res$coe[a])
    rcoe = c(rcoe, res$coe[paste0("sd", sep = ".", a)])
  } else if (res$ranp[a] == "n") {
    typepar[length(kcoe) + 1] = "n"
    coe = c(coe, res$coe[a])
    rcoe = c(rcoe, res$coe[paste0("sd", sep = ".", a)])
  }
  names(typepar) = c(names(res$coe)[c(kcoe, costcoe)])
  coe = c(coe, rcoe)

  # Draws values of the random parameters from their appropriate distributions.
  # Creates a list beta with a number of elements equal to the number of parameters.
  # Each element containes the Q draws from the distribution of that parameter
  # (if random) or the unique value of the parameter (if fixed).
  # Also register the random draws from normal
  beta = list()
  Z = list()
  for (i in (1:ncoe))
  {
    if (typepar[i] == "fixed") {
      beta = c(beta, list(res$coe[names(typepar[i])]))
      Z = c(Z, 0)
    } else if (typepar[i] == "ln") {
      z = rnorm(Q)
      beta = c(beta, list(exp(coe[names(typepar[i])] + coe[paste0("sd", sep =
                                                                    ".", names(typepar[i]))] * z)))
      Z = c(Z, list(z))
    } else if (typepar[i] == "n") {
      z = rnorm(Q)
      beta = c(beta, list(coe[names(typepar[i])] + coe[paste0("sd", sep =
                                                                ".", names(typepar[i]))] * z))
      Z = c(Z, list(z))
    }
  }

  # If the cost parameter follows a lognormal distribution,
  # the sign of beta_c is corrected
  # so that WTP has the right sign in the output
  if (typepar[ncoe] == "ln")
    beta[[ncoe]] = -beta[[ncoe]]

  # Compute the conditional WTP (equation 9 in the paper)
  mDWTP_q = matrix(0, Q, length(kcoe))
  for (i in (1:length(kcoe)))
  {
    mDWTP_q[, i] = -beta[[i]] / beta[[ncoe]]
  }

  ## Compute point estimate of the average WTP, using a function of the
  # structural parameter or the average of the conditional WTP
  # (equation 12 in the paper). In case the cost parameter follows a gaussian
  # distribution, the median WTP is provided instead of the mean.
  mDWTP = rep(0, ncoe - 1)
  for (k in 1:(ncoe - 1))
  {
    if (typepar[k] == "n" & typepar[ncoe] == "fixed") {
      mDWTP[k] = -coe[k] / coe[ncoe]
    } else if (typepar[k] == "fixed" & typepar[ncoe] == "ln") {
      mDWTP[k] = sign(coe[k]) * exp(-coe[ncoe] + log(abs(coe[k])) + coe[length(coe)] ^
                                      2 / 2)
    } else if (typepar[k] == "ln" & typepar[ncoe] == "ln") {
      mDWTP[k] = exp(coe[k] - coe[ncoe] + (rcoe[k] ^ 2 + rcoe[ncoe] ^ 2) / 2)
    } else {
      mDWTP[k] = mean(mDWTP_q[, k])
    }
  }
  if (typepar[ncoe] == "n") {
    mDWTP = apply(mDWTP_q, 2, median)
  }
  names(mDWTP) = names(res$coe)[kcoe]

  # Computes standard error of conditional WTP (as in equation 14 of the paper)
  # Then computes standard error of the mean WTP (as in equation 15 of the paper)
  covDWTP_q = array(0, dim = c(Q, Q, ncoe - 1))
  sebisDWTP_q = matrix(0, Q, length(kcoe))
  semean = rep(0, length(kcoe))

  for (i in 1:length(kcoe))
  {
    if (typepar[i] == "fixed")
    {
      if (typepar[ncoe] == "ln")
      {
        fdbis = rbind(-1 / beta[[ncoe]], beta[[i]] / beta[[ncoe]], beta[[i]] / beta[[ncoe]] *
                        Z[[ncoe]])
        covDWTP_q[, , i] = t(fdbis) %*% (vcov(res)[c(names(res$coe)[kcoe[i]],
                                                     names(res$coe)[costcoe],
                                                     paste0("sd", sep = '.', names(res$coe)[costcoe])), c(names(res$coe)[kcoe[i]],
                                                                                                          names(res$coe)[costcoe],
                                                                                                          paste0("sd", sep =
                                                                                                                   '.', names(res$coe)[costcoe]))]) %*% fdbis
        sebisDWTP_q[, i] = sqrt(diag(covDWTP_q[, , i]))
      }
    } else if (typepar[i] == "ln")
    {
      if (typepar[ncoe] == "ln")
      {
        fdbis = rbind(-beta[[i]] / beta[[ncoe]],
                      beta[[i]] / beta[[ncoe]],-Z[[i]] * beta[[i]] / beta[[ncoe]],
                      Z[[ncoe]] * beta[[i]] / beta[[ncoe]])
        covDWTP_q[, , i] = t(fdbis) %*% vcov(res)[c(
          names(res$coe)[kcoe[i]],
          names(res$coe)[costcoe],
          paste0("sd", sep = '.', names(res$coe)[kcoe[i]]),
          paste0("sd", sep = '.', names(res$coe)[costcoe])
        ), c(
          names(res$coe)[kcoe[i]],
          names(res$coe)[costcoe],
          paste0("sd", sep =
                   '.', names(res$coe)[kcoe[i]]),
          paste0("sd", sep = '.', names(res$coe)[costcoe])
        )] %*% fdbis
        sebisDWTP_q[, i] = sqrt(diag(covDWTP_q[, , i]))
      }
    } else if (typepar[i] == "n")
    {
      if (typepar[ncoe] == "fixed")
      {
        fdbis = rbind(-1 / beta[[ncoe]], beta[[i]] / beta[[ncoe]] ^ 2,-Z[[i]] /
                        beta[[ncoe]])
        covDWTP_q[, , i] = t(fdbis) %*% vcov(res)[c(names(res$coe)[kcoe[i]],
                                                    names(res$coe)[costcoe],
                                                    paste0("sd", sep = '.', names(res$coe)[kcoe[i]])), c(names(res$coe)[kcoe[i]],
                                                                                                         names(res$coe)[costcoe],
                                                                                                         paste0("sd", sep =
                                                                                                                  '.', names(res$coe)[kcoe[i]]))] %*% fdbis
        sebisDWTP_q[, i] = sqrt(diag(covDWTP_q[, , i]))
      } else if (typepar[ncoe] == "n") {
        fdbis = rbind(-1 / beta[[ncoe]], beta[[i]] / beta[[ncoe]] ^ 2,-Z[[i]] /
                        beta[[ncoe]], beta[[i]] / beta[[ncoe]] ^ 2 * Z[[ncoe]])
        covDWTP_q[, , i] = t(fdbis) %*% vcov(res)[c(
          names(res$coe)[kcoe[i]],
          names(res$coe)[costcoe],
          paste0("sd", sep = '.', names(res$coe)[kcoe[i]]),
          paste0("sd", sep = '.', names(res$coe)[costcoe])
        ), c(
          names(res$coe)[kcoe[i]],
          names(res$coe)[costcoe],
          paste0("sd", sep =
                   '.', names(res$coe)[kcoe[i]]),
          paste0("sd", sep = '.', names(res$coe)[costcoe])
        )] %*% fdbis
        sebisDWTP_q[, i] = sqrt(diag(covDWTP_q[, , i]))
      }
    }

    semean[i] = sqrt(mean(covDWTP_q[, , i]))
  }

  # Computes prediction error (as in equations 17, 18 and 19 in the paper)
  setot = sqrt(colMeans(sebisDWTP_q ^ 2) + apply(mDWTP_q, 2, var))

  # Computes confidence interval for the mean of WTP (using normal theory and
  # equation 15 in the paper).
  CIm = list()
  for (i in (1:length(level)))
  {
    CIm = c(CIm, list(rbind(
      mDWTP + qnorm((1 - level[i]) / 2) * semean, mDWTP + qnorm(level[i] + (1 -
                                                                              level[i]) / 2) * semean
    )))
    names(CIm)[i] = paste0("confidence level", sep = "=", level[i])
    rownames(CIm[[i]]) = c("lower bound", "upper bound")
  }

  # Computes prediction interval for WTP (using equations 20 and 21 in the paper)
  CIbis = list()
  cik = matrix(0, nrow = 2, ncol = (ncoe - 1))
  for (i in (1:length(level)))
  {
    for (k in (1:(ncoe - 1)))
    {
      cik[, k] = c(
        F_inv((1 - level[i]) / 2,
              1 / Q,
              mDWTP_q[, k],
              sebisDWTP_q[, k],
              br = c(-1000, 1000)
        ),
        F_inv(
          1 - (1 - level[i]) / 2,
          1 / Q,
          mDWTP_q[, k],
          sebisDWTP_q[, k],
          br = c(-1000, 1000)
        )
      )
    }
    CIbis = c(CIbis, list(cik))
  }


  ## estimated tail probabilities
  tDWTP = matrix(rep(0, (ncoe - 1)), 1, ncoe - 1)
  for (j in (1:(ncoe - 1)))
  {
    tDWTP[j] = mean(mDWTP_q[, j] > qw)
  }


  ## estimated quantiles
  qDWTP = apply(mDWTP_q, 2, quantile, probs = 1 - tpw, na.rm = TRUE)

  results = list(mDWTP, semean, setot, CIm, CIbis, qDWTP, tDWTP, wtp_names)
  names(results) = c(
    "WTP",
    "se.mean.WTP",
    "pred.err.WTP",
    "CI.mean.WTP",
    "PI.WTP",
    "quant.WTP",
    "tail.prob.WTP",
    "wtpnames"
  )
  return(results)
}


#############Krinsky & Robb method###################################

#' Krinsky & Robb
#'
#' @param res an object of class gmnl
#' @param Q number of random draws from standard distributions (normal or uniform)
#' @param B number of bootstrap samples
#' @param kcoe vector of positions (in the res$coe vector) of the parameters for which wtp is to be computed
#' @param costcoe position (in the res$coe vector) of the cost parameter
#' @param level vector of confidence level to compute confidence and prediction intervals (i.e. 0.9, 0.95…)
#' @param qw value (the same for all parameter for which wtp is to be computed) at which upper tail probability
#'           for wtp distribution is to be computed
#' @param tpw probability (the same for all parameters for which wtp is to be computed) at which quantiles leaving
#'            that probability in the upper tail of wtp distribution are to be computed
#'
#' @importFrom MASS mvrnorm
#'
#' @return vector of : "se.mean.WTP","pred.err.WTP","CI.mean.WTP","PI.WTP","quant.WTP","se.quant.WTP","CI.quant.WTP","tail.prob.WTP","se.tail.prob.WTP","CI.tail.prob.WTP"
#' @export

KRCI <- function(res, Q, B, kcoe, costcoe, level, qw, tpw) {
  ncoe = length(kcoe) + length(costcoe)
  coe = c()
  rcoe = c()
  wtp_names <- res$wtp_names

  # Correcting for negative estimates of standard deviation
  # for random coefficients (sometimes happening with mlogit,
  # not sure this is required for gmnl)
  rp = length(res$ranp)
  for (i in (1:rp))
  {
    res$coe[paste0("sd", sep = '.', names(res$ranp)[i])] = abs(res$coe[paste0("sd", sep =
                                                                                '.', names(res$ranp)[i])])
  }

  # Determins the type of distribution for paramaters of interest
  # and creates an ordered vector with coefficients (or means if random) of the attributes for which
  # wtp needs to be calculated, the cost coefficient (or mean if random)
  # and standard deviations of random coefficients
  typepar = rep(0, ncoe)
  for (i in (1:length(kcoe)))
  {
    a = names(res$coefficients)[kcoe[i]]
    if (sum(a == names(res$ranp)) == 0) {
      typepar[i] = "fixed"
      coe = c(coe, res$coefficients[a])
    } else if (res$ranp[a] == "ln") {
      typepar[i] = "ln"
      coe = c(coe, res$coefficients[a])
      rcoe = c(rcoe, res$coefficients[paste0("sd", sep = ".", a)])
    } else if (res$ranp[a] == "n") {
      typepar[i] = "n"
      coe = c(coe, res$coefficients[a])
      rcoe = c(rcoe, res$coefficients[paste0("sd", sep = ".", a)])
    }
  }
  a = names(res$coefficients)[costcoe]
  if (sum(a == names(res$ranp)) == 0) {
    typepar[length(kcoe) + 1] = "fixed"
    coe = c(coe, res$coefficients[a])
  } else if (res$ranp[a] == "ln") {
    typepar[length(kcoe) + 1] = "ln"
    coe = c(coe, res$coefficients[a])
    rcoe = c(rcoe, res$coefficients[paste0("sd", sep = ".", a)])
  } else if (res$ranp[a] == "n") {
    typepar[length(kcoe) + 1] = "n"
    coe = c(coe, res$coefficients[a])
    rcoe = c(rcoe, res$coefficients[paste0("sd", sep = ".", a)])
  }
  names(typepar) = c(names(res$coefficients)[c(kcoe, costcoe)])
  coe = c(coe, rcoe)

  # First level draws from multivariate normal to account for
  # estimation error (draws of the structural parameters or the fixed parameters)
  KRlevel1 <-
    MASS::mvrnorm(
      B,
      coe,
      vcov(res)[names(coe), names(coe)], # this was formerly vcov(res)[names(coe), names(coe)]
      tol = 1e-6,
      empirical = FALSE,
      EISPACK = FALSE
    )

  # Second level draws from the distribution of each coefficient from the respective mixing distribution
  KRlevel2 <- array(0, dim = c(B, Q, ncoe))
  pos = which(typepar != "fixed")
  z = rnorm(Q * rp)
  z = matrix(z, Q, rp)
  for (i in 1:B)
  {
    for (k in 1:ncoe)
    {
      if (typepar[k] == "fixed") {
        KRlevel2[i, , k] = KRlevel1[i, k]
      } else if (typepar[k] == "ln") {
        KRlevel2[i, , k] = exp(KRlevel1[i, k] + KRlevel1[i, paste0("sd", sep = ".", names(coe[k]))] *
                                 z[, which(pos == k)])
      } else if (typepar[k] == "n") {
        KRlevel2[i, , k] = KRlevel1[i, k] + KRlevel1[i, paste0("sd", sep = ".", names(coe[k]))] *
          z[, which(pos == k)]
      }
    }
  }

  # If the cost parameter is ln, the sign of beta_c is corrected so that WTP has the right sign in the output
  if (typepar[ncoe] == "ln")
    KRlevel2[, , k] = -KRlevel2[, , k]

  # Compute the Q*B*(ncoe-1) simulated wpt
  # Results are in an array of dimension BxQxkcoe so that for b=1 and kcoe = 1
  # we have the distribution of KWTP1 due to heterogeneity while for q=1 and kcoe=1
  # we have the distribution of KWTP1 due to standard error
  KWTP = -KRlevel2[, , 1:(ncoe - 1)] / array(KRlevel2[, , ncoe], dim = c(B, Q, length(kcoe)))

  ## point estimate of the average WTP
  # and of the average WTP in each bootstrap sample,
  # to give bootstrap standard errors for the mean WTP
  # In case the cost parameter follows a gaussian
  # distribution, the median WTP is provided instead of the mean.
  mKWTP = rep(0, ncoe - 1)
  mKWTP_b = matrix(0, B, ncoe - 1)
  for (k in 1:(ncoe - 1)) {
    if (typepar[k] == "n" & typepar[ncoe] == "fixed") {
      mKWTP[k] = -coe[k] / coe[ncoe]
      mKWTP_b[, k] = -KRlevel1[, k] / KRlevel1[, ncoe]
    } else if (typepar[k] == "fixed" & typepar[ncoe] == "ln") {
      mKWTP[k] = sign(coe[k]) * exp(-coe[ncoe] + log(abs(coe[k])) + coe[length(coe)] ^
                                      2 / 2)
      mKWTP_b[, k] = sign(KRlevel1[, k]) * exp(-KRlevel1[, ncoe] + log(abs(KRlevel1[, k])) +
                                                 KRlevel1[, ncoe + 1] ^ 2 / 2)
    } else if (typepar[k] == "ln" & typepar[ncoe] == "ln") {
      mKWTP[k] = exp(coe[k] - coe[ncoe] + (coe[ncoe + k] ^ 2 + coe[length(coe)] ^
                                             2) / 2)

      mKWTP_b[, k] = exp(KRlevel1[, k] - KRlevel1[, ncoe] + (KRlevel1[, ncoe +
                                                                        k] ^ 2 + KRlevel1[, length(coe)] ^ 2) / 2)
    } else {
      mKWTP[k] = mean(KWTP[, , k])
      mKWTP_b[, k] =	apply(KWTP[, , k], 1, mean)
    }
  }
  if (typepar[ncoe] == "n") {
    mKWTP = apply(KWTP, 3, median)
  }

  # Standard error of the average WTP (due to sampling error),
  # computed as standard deviations of bootstrap mean WTPs
  semediecond_b = apply(mKWTP_b, 2, sd)

  # Confidence interval for the mean
  CImKWTP = apply(mKWTP_b,
                  2,
                  quantile,
                  probs = c((1 - level) / 2, level + (1 - level) / 2),
                  na.rm = TRUE)

  # Estimated prediction error
  peKWTP = apply(KWTP, 3, sd)

  ## prediction interval for WTP
  CI = apply(KWTP,
             3,
             quantile,
             probs = c((1 - level) / 2, level + (1 - level) / 2),
             na.rm = TRUE)

  ## estimated tail probabilities
  tKWTP_b = matrix(rep(0, (ncoe - 1) * B), B, ncoe - 1)
  for (j in (1:(ncoe - 1)))
    tKWTP_b[, j] = rowMeans(KWTP[, , j] > qw)
  tp = colMeans(tKWTP_b)

  ## estimated standard error for tail probabilities
  setp = apply(tKWTP_b, 2, sd)

  ## Confidence interval for tail probabilities
  CItp = apply(tKWTP_b,
               2,
               quantile,
               probs = c((1 - level) / 2, level + (1 - level) / 2),
               na.rm = TRUE)

  ## estimated quantiles
  qKWTP = apply(KWTP, 3, quantile, probs = 1 - tpw, na.rm = TRUE)

  ## estimated standard error for estimated quantiles
  qKWTP_b = apply(KWTP,
                  c(1, 3),
                  quantile,
                  probs = 1 - tpw,
                  na.rm = TRUE)
  seQ = apply(qKWTP_b, 2, sd)

  ## CI for estimated quantiles
  CIq = apply(qKWTP_b,
              2,
              quantile,
              probs = c((1 - level) / 2, level + (1 - level) / 2),
              na.rm = TRUE)

  results = list(mKWTP,
                 semediecond_b,
                 peKWTP,
                 CImKWTP,
                 CI,
                 qKWTP,
                 seQ,
                 CIq,
                 tp,
                 setp,
                 CItp,
                 wtp_names)
  names(results) = c(
    "WTP",
    "se.mean.WTP",
    "pred.err.WTP",
    "CI.mean.WTP",
    "PI.WTP",
    "quant.WTP",
    "se.quant.WTP",
    "CI.quant.WTP",
    "tail.prob.WTP",
    "se.tail.prob.WTP",
    "CI.tail.prob.WTP",
    "wtp_names"
  )
  return(results)
}


#' Summary function for DeltaCI Method Output
#'
#' @param object output object from DeltaCI function
#' @param ... other arguments
#'
#' @return visual table of the results
#' @export
summary.Delta <- function(object, ...) {
  wtp_names <- object$wtpnames
  n <- length(wtp_names)

  # Define readable metric names
  metric_names <- c("Mean WTP", "SE Mean WTP", "Pred. Err. WTP", "CI Mean WTP",
                    "PI WTP", "Quant. WTP","Tail Prob. WTP")

  # Print header row with coefficient names
  cat(strrep(" ", 17), " | ")
  cat(paste(sprintf("%-17s", wtp_names), collapse = " | "), "\n")

  # Print separator
  cat(strrep("-", 20 + 20 * n), "\n")

  # Print metric names and corresponding values
  for (i in 1:length(metric_names)) {
    cat(sprintf("%-18s | ", metric_names[i]))

    # Check if the current metric is numeric or character
    if (is.numeric(unlist(object[i])[1])) {
      if (length(unlist(object[i])) %% 2 == 0 && i %in% c(4, 5, 8, 11)) {
        formatted_intervals <- sapply(seq(1, length(unlist(object[i])), 2), function(j) {
          sprintf("[%7.3f,%7.3f]", unlist(object[i])[j], unlist(object[i])[j + 1])
        })
        cat(paste(formatted_intervals, collapse = " | "), "\n")
      } else {
        cat(paste(sprintf("% 17.8f", unlist(object[i])[1:n]), collapse = " | "), "\n")
      }
    } else {
      cat(paste(sprintf("%17s", unlist(object[i])[1:n]), collapse = " | "), "\n")
    }
  }
}


#' Summary function for KRCI
#'
#' @param object KRCI object
#' @param ... other
#'
#' @export
summary.KRCI <- function(object, ...) {
  wtp_names <- object$wtp_names
  n <- length(wtp_names)

  # Define readable metric names
  metric_names <- c("Mean WTP", "SE Mean WTP", "Pred. Err. WTP", "CI Mean WTP",
                    "PI WTP", "Quant. WTP", "SE Quant. WTP", "CI Quant. WTP",
                    "Tail Prob. WTP", "SE Tail Prob. WTP", "CI Tail Prob. WTP")

  # Print header row with coefficient names
  cat(strrep(" ", 17), " | ")
  cat(paste(sprintf("%-17s", wtp_names), collapse = " | "), "\n")

  # Print separator
  cat(strrep("-", 20 + 20 * n), "\n")

  # Print metric names and corresponding values
  for (i in 1:length(metric_names)) {
    cat(sprintf("%-18s | ", metric_names[i]))

    # Check if the current metric is numeric or character
    if (is.numeric(unlist(object[i])[1])) {
      if (length(unlist(object[i])) %% 2 == 0 && i %in% c(4, 5, 8, 11)) {
        formatted_intervals <- sapply(seq(1, length(unlist(object[i])), 2), function(j) {
          sprintf("[%7.3f,%7.3f]", unlist(object[i])[j], unlist(object[i])[j + 1])
        })
        cat(paste(formatted_intervals, collapse = " | "), "\n")
      } else {
        cat(paste(sprintf("% 17.8f", unlist(object[i])[1:n]), collapse = " | "), "\n")
      }
    } else {
      cat(paste(sprintf("%17s", unlist(object[i])[1:n]), collapse = " | "), "\n")
    }
  }
}

# Define custom estfun function for gmnl models
#'
#' @param model of gmnl model class
#'
#' @return matrix
estfun.gmnl <- function(model) {
  return(model$logLik$gradientObs)
}

#' Compute Clustered Standard Errors for a gmnl Model
#'
#' Computes the clustered standard errors for a gmnl model object
#' using the sandwich `vcovCL()` function.
#'
#' @param model   A `gmnl` model object for which the clustered standard errors will be computed.
#' @return        Prints the clustered standard errors covariance matrix.
#' @importFrom sandwich vcovCL
#' @export
compute_clustered_standard_error.gmnl <- function(model) {
  # if (!inherits(model, "gmnl")) {
  if (toString(class(model)) != "gmnl") {
    stop("The input model must be of class 'gmnl'")
  }

  print(sandwich::vcovCL(model))
}

#' Summation of Log of Inverse Number of Alternatives for Each Observation
#'
#' This function calculates the summation of the log of the inverse of the number of alternatives
#' for each observation in the data. Used in calculating the null model for McFadden R^2.
#'
#' @param model A gmnl model object.
#' @return A numeric value representing the summation of the log of the inverse of the number
#'   of alternatives for each observation.
sum_log_inv_num_alternatives <- function(model) {
  if (!inherits(model, "gmnl")) {
    stop("The input model must be of class 'gmnl'")
  }

  # Extract the number of alternatives for each observation from the data
  num_alternatives <- tapply(model$data$alt, INDEX = model$data$chid, FUN = length)

  # Calculate the log of the inverse of the number of alternatives for each observation
  log_inv_num_alternatives <- log(1 / num_alternatives)

  # Sum the result
  sum_log_inv_num_alternatives <- sum(log_inv_num_alternatives)

  return(sum_log_inv_num_alternatives)
}


#' Compute elasticities for a GMNL model
#'
#' Given a GMNL model and a reference alternative, this function computes the
#' elasticity of each attribute with respect to the reference alternative. The
#' function returns a summary of the elasticities by mode.
#'
#' @param gmnl_model A GMNL model object
#' @param wrta A character string specifying the reference alternative with
#' respect to which the elasticities are to be computed.
#' @param cost Name of cost variable
#'
#' @return A gmnl object with the updated elasticities, along with the data summary
#'
#' @importFrom stats model.matrix fitted.values
#' @importFrom dplyr group_by summarise %>% .data
#' @export
compute_elasticities_direct.gmnl <- function(gmnl_model, wrta) {
  # if (!inherits(gmnl_model, "gmnl")) {
  if (toString(class(gmnl_model)) != "gmnl") {
    stop("The input model must be of class 'gmnl'")
  }

  # Get the data from the gmnl_model object
  data <- gmnl_model$data

  # Get the number of alternatives
  num_alts <- length(gmnl_model$freq)

  # Predict using inbuilt function
  predicted_probabilities <- fitted.values(gmnl_model, outcome = FALSE)
  data$pred_prob <- as.vector(t(predicted_probabilities))

  # Predictors used in the model
  X <- model.matrix(gmnl_model)
  exp_product <- as.numeric(exp(X %*% coef(gmnl_model)))

  # k-th predictor: 1st predictor (the changing attribute)
  x_kij <- as.numeric(X[, "vcost"])

  # Creating (N x K) matrices, where K is the number of alternatives
  x_matrix <- matrix(x_kij, ncol = num_alts, byrow = TRUE)
  exp_product_matrix <- matrix(exp_product, ncol = num_alts, byrow = TRUE)

  # k-th alternative
  wrta_index <- which(colnames(gmnl_model$prob.alt) == wrta)
  exp_product_l <- exp_product_matrix[, wrta_index]
  x_kil <- x_matrix[, wrta_index]

  exp_product_matrix_sum <- apply(exp_product_matrix, 1, sum)
  exp_product_matrix_sum <- rep(exp_product_matrix_sum, each = num_alts)
  exp_product_l <- rep(exp_product_l, each = num_alts)
  x_kil <- rep(x_kil, each = num_alts)

  # Compute the choice probabilities
  p_ij <- exp_product / exp_product_matrix_sum
  p_il <- exp_product_l / exp_product_matrix_sum

  # Compute the elasticities
  data$elast1[data$mode == wrta] <- (1 - p_ij[data$mode == wrta]) * coef(gmnl_model)["vcost"] * x_kij[data$mode == wrta]
  data$elast1[data$mode != wrta] <- (-p_il[data$mode != wrta]) * coef(gmnl_model)["vcost"] * x_kil[data$mode != wrta]

  # Use dplyr to summarize the data by mode
  data_summary <- data %>%
    group_by(mode) %>%
    summarise(actual = mean(choice),
              predicted = mean(pred_prob),
              elasticity = mean(elast1))
  # Update the data in the gmnl_model object with the computed elasticities
  gmnl_model$data <- data

  # Return the updated gmnl_model object with the computed elasticities
  return(list(model = gmnl_model, data_summary = data_summary))
}


#' Compute elasticities for a given generalized multinomial logit (gmnl) model
#'
#' This function computes the cross elasticities for a given gmnl model by changing the value of the cost variable for a specific alternative and then computing the change in predicted probabilities.
#'
#' @param gmnl_model A gmnl model object.
#' @param wrta A string indicating the alternative to compute elasticities for.
#' @param price Name of the cost variable
#'
#' @return A gmnl model object with updated elasticities and predicted probabilities.
#'
#' @importFrom stats model.matrix fitted.values
#' @importFrom dplyr group_by summarise %>% .data
#' @export
compute_elasticities_cross.gmnl <- function(gmnl_model, wrta) {
  # if (!inherits(gmnl_model, "gmnl")) {
  if (toString(class(gmnl_model)) != "gmnl") {
    stop("The input model must be of class 'gmnl'")
  }

  # Get the data from the gmnl_model object
  data <- gmnl_model$data

  # Get the number of alternatives
  num_alts <- length(gmnl_model$freq)

  # Update the cost variable for the specified alternative
  data$vcost[data$mode == wrta] <- data$vcost[data$mode == wrta] * 1.01

  # Predict using inbuilt function
  Xnew <- model.matrix(gmnl_model)
  Xnew[,"vcost"] <- as.matrix(data$vcost)
  eXb <- as.numeric(exp(Xnew %*% coef(gmnl_model)))
  SeXb <-  matrix(eXb, ncol = num_alts, byrow = TRUE)
  SeXb <-  apply(SeXb, 1, sum)
  SeXb <-  rep(SeXb, each = num_alts)
  pr <- eXb / SeXb
  data$pred_prob_new <- pr

  # Compute the change in predicted probabilities and elasticities
  data$prob_chg <- 100 * (data$pred_prob_new / data$pred_prob - 1)
  data$elast2 <- data$prob_chg / 1

  data_summary <- data %>%
    group_by(mode) %>%
    summarise(actual = mean(choice),
              predicted = mean(pred_prob),
              change = mean(prob_chg),
              direct_elast = mean(elast1),
              cross_elast = mean(elast2))

  # Update the elasticities in the gmnl_model object
  gmnl_model$data <- data

  return(list(model = gmnl_model, data_summary = data_summary))
}

