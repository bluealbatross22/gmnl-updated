# Compute time
make.time <- function(object){
  et <- object$time[3]
  s <- round(et,0)
  h <- s%/%3600 # hours
  s <- s-3600*h
  m <- s%/%60 # minutes
  s <- s-60*m
  paste(h, "h:", m, "m:", s, "s", sep = "")
}

# Returns n columns of matrix/array turned into row vector(s)
repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

# Returns n columns of matrix/array turned into column vector(s)
repCols <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

# makeL() is own function to make lower triangular

# Constructs a polynomial term for a specific row and column in an additive
# covariance matrix.
#
# @param row The row index in the covariance matrix.
# @param col The column index in the covariance matrix.
# @param Ka The dimension of the covariance matrix.
#
# @return A string representation of the polynomial term for the specified row
#         and column in the additive covariance matrix.
make.add <- function(row, col, Ka) {
  # Validate input parameters
  if (!is.numeric(row) || !is.numeric(col) || !is.numeric(Ka)) {
    stop("row, col, and Ka must be numeric.")
  }

  if (row < 1 || col < 1 || Ka < 1) {
    stop("row, col, and Ka must be greater than or equal to 1.")
  }

  # Create an nxn lower triangular matrix using makeL function (assumes makeL is defined elsewhere)
  lower_triangular_matrix <- makeL(1:rep(0.5 * Ka * (Ka + 1)))

  polynomial_term <- ""
  for (k in row:col) {
    curr_base_index <- lower_triangular_matrix[k, row]

    # Update the polynomial term by concatenating variable names with indices
    indices <- curr_base_index:(curr_base_index + (Ka - k))
    polynomial_term <- paste(paste("x", indices, sep = ""), paste("x", curr_base_index, sep = ""), sep = "*")
  }

  # Return the polynomial term as a string
  return(polynomial_term)
}


# Returns an array/matrix containing n copies of A in the row and column dimensions
# If x is a matrix it appends replicas to the bottom
# If x is a vector it appends to the right
repmat <- function(x, dimen){
  kronecker(array(1, dim = dimen), x)
}


