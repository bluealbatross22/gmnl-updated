#' GMNL Data
#'
#' @importFrom dfidx dfidx
#' @param data a `data.frame`,
#' @param choice the variable indicating the choice made: it can be
#'     either a logical vector, a numerical vector with 0 where the
#'     alternative is not chosen, a factor with level 'yes' when the
#'     alternative is chosen
#' @param shape the shape of the `data.frame`: whether `long` if each
#'     row is an alternative or `wide` if each row is an observation
#' @param varying the indexes of the variables that are alternative
#'     specific, sometimes seen as xxxx.1, xxxx.2, etc
#' @param sep the charachter  that separates the variable name
#'     and the alternative name (only for `wide` data)
#' @param alt.var the name of the variable that contains the
#'     alternative index (for a `long` `data.frame` only) or the name
#'     under which the alternative index will be stored (the default
#'     name is `alt`),
#' @param chid.var the name of the variable that contains the choice
#'     index or the name under which the choice index will be stored,
#' @param alt.levels the name of the alternatives: if null, for a
#'     `wide` data.frame, they are guessed from the variable names and
#'     the choice variable (both should be the same), for a `long`
#'     `data.frame`, they are guessed from the `alt.var` argument,
#' @param id.var the name of the variable that contains the individual
#'     index if any,
#' @param group.var the name of the variable that contains the group
#'     index if any,
#' @param opposite returns the opposite of the specified variables,
#' @param drop.index should the index variables be dropped from the
#'     `data.frame`
#' @param ranked a logical value which is true if the response is a ranked variable
#' @param subset logical expression which defines the subset of
#'     observations to be selected
#' @param ... more arguments to specify reshaping of the data
#'
#' @return gmnl.data object that inherits class dfidx
#' @export
#'
gmnl.data <- function (data, choice = NULL, shape = c("long", "wide"), varying = NULL,
                         sep = ".", alt.var = NULL, chid.var = NULL, alt.levels = NULL,
                         id.var = NULL, group.var = NULL, opposite = NULL, drop.index = FALSE,
                         ranked = FALSE, subset = NULL, ...) {
  call_data <- match.call(expand.dots = TRUE)

  # Handle the 'subset' parameter
  data <- handle_subset(data, call_data)
  idx <- setup_idx(alt.var, chid.var, group.var, id.var)

  call_data$idx <- idx
  call_data$idnames <- c("chid", "alt")
  call_data$drop.index <- drop.index

  # Rename 'alt.levels' to 'levels' in call_data
  alt_levels_pos <- match("alt.levels", names(call_data))
  if (!is.na(alt_levels_pos)) names(call_data)[alt_levels_pos] <- "levels"

  # Extract relevant names from call_data
  relevant_names <- c("data", "idx", "choice", "shape", "varying", "sep", "levels", "opposite", "idnames",
                      "pkg", "drop.index", "subset", "ranked")
  matched_positions <- match(relevant_names, names(call_data), 0)
  call_data <- call_data[c(1, matched_positions)]
  call_data[[1]] <- as.name("dfidx")

  # Set up default argument values for dfidx function
  dfidx_defaults <- list(data = NA, idx = NULL, drop.index = TRUE, as.factor = NULL, pkg = NULL,
                         fancy.row.names = FALSE, subset = NULL,
                         idnames = NULL, shape = "long", choice = NULL,
                         varying = NULL, sep = ".", opposite = NULL, levels = NULL, ranked = FALSE)

  # Replace default values with those specified in gmnl.data
  dfidx_args <- dfidx_defaults
  dfidx_args[names(call_data)[-1]] <- as.list(call_data)[-1]

  # Evaluate names in dfidx_args if needed
  for (i in 2:length(dfidx_args)) {
    if (is.name(dfidx_args[[i]])) {
      dfidx_args[[i]] <- eval(dfidx_args[[i]], parent.frame())
    }
  }

  gmnl_data <- dfidx::dfidx(data = data, dfidx_args$idx, drop.index = dfidx_args$drop.index,
                            as.factor = dfidx_args$as.factor,
                            fancy.row.names = dfidx_args$fancy.row.names,
                            idnames = dfidx_args$idnames, shape = dfidx_args$shape,
                            choice = dfidx_args$choice, varying = dfidx_args$varying,
                            sep = dfidx_args$sep, opposite = dfidx_args$opposite,
                            levels = dfidx_args$levels, ranked = dfidx_args$ranked)

  # Add 'gmnl.data' as an additional class to the returned object
  class(gmnl_data) <- c(class(gmnl_data), "gmnl.data")
  return(gmnl_data)
}

# Helper function to handle the 'subset' parameter
handle_subset <- function(data, cldata) {
  if ("subset" %in% names(cldata)) {
    sub_call <- cldata
    m <- match(c("data", "subset"), names(sub_call), 0)
    sub_call <- sub_call[c(1, m)]
    names(sub_call)[2] <- "x"
    sub_call[[1]] <- as.name("subset")
    data <- eval(sub_call, parent.frame())
  }
  return(data)
}

# Helper function to set up the 'idx' variable
setup_idx <- function(alt_var, chid_var, group_var, id_var) {
  idx <- NULL

  if (is.null(alt_var) & is.null(chid_var)) idx <- NULL
  if (is.null(alt_var) & !is.null(chid_var)) idx <- list(chid_var)
  if (!is.null(alt_var) & is.null(chid_var)) idx <- list(NA, alt_var)
  if (!is.null(alt_var) & !is.null(chid_var)) idx <- list(chid_var, alt_var)

  if (!is.null(group_var)) idx[[2]] <- c(idx[[2]], group_var)
  if (!is.null(id_var)) {
    if (is.null(idx)) idx <- list(c(NA, id_var))
    else idx[[1]] <- c(idx[[1]], id_var)
  }

  return(idx)
}


#gFormula from current version of mlogit
model.matrix.gFormula2 <- function(object, data, ...){
  if (inherits(data, "dfidx")){
    mf <- model.frame(formula = data, data = object, ...)
    stats::model.matrix(mf, ...)
  }
  else{
    class(object) <- c("Formula", "formula")
    model.matrix(object, data, ...)
  }
}

# Time series object library
#' @import zoo
zoo::index

# Indexed data frame data type library
#' @importFrom dfidx idx
index.dfidx <- function(x, ...){
  dfidx::idx(x, ...)
}

#' @importFrom dfidx idx
index.gmnl <- function(x, ...){
  dfidx::idx(x, ...)
}

