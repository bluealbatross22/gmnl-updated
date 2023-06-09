% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmnl.extras.R
\name{sum_log_inv_num_alternatives}
\alias{sum_log_inv_num_alternatives}
\title{Summation of Log of Inverse Number of Alternatives for Each Observation}
\usage{
sum_log_inv_num_alternatives(model)
}
\arguments{
\item{model}{A gmnl model object.}
}
\value{
A numeric value representing the summation of the log of the inverse of the number
  of alternatives for each observation.
}
\description{
This function calculates the summation of the log of the inverse of the number of alternatives
for each observation in the data. Used in calculating the null model for McFadden R^2.
}
