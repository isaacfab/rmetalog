#' rmetalog: R implementation of the metalog distribution
#'
#' The \code{rmetalog} package implements the metalog distribution in \code{R}
#'
#' @docType package
#' @name rmetalog
#' @importFrom ggplot2 aes

if(getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      "cumValue",
      "pdfValues",
      "quantileValues"
    )
  )


#' Fit the metalog distribution to data
#'
#' @param x vector of numeric data
#' @param term_limit integer between 3 and 30, specifying the number of metalog
#'   distributions to generate. Larger term distributions have more flexibility
#'   (default: 13)
#' @param bounds numeric vector specifying lower or upper bounds, none required
#'   if the distribution is unbounded
#' @param boundedness character string specifying unbounded, semi-bounded upper,
#'   semi-bounded lower or bounded; accepts values \code{u}, \code{su},
#'   \code{sl} and \code{b} (default: 'u')
#' @param term_lower_bound (Optional) the smallest term to generate, used to
#'   minimize computation of unwanted terms must be less than term_limit (default is 2)
#' @param step_len (Optional) size of steps to summarize the distribution
#'   (between 0 and 0.01) this is only used if the data vector length is greater
#'   than 100. Use this if a specific fine grid fit is required. (default is
#'   0.01)
#' @param probs (Optional) probability quantiles, same length as \code{x}
#' @param fit_method (Optional) preferred method of fitting distribution: accepts values
#'  \code{OLS}, \code{LP} or \code{any} (defaults to any)
#'
#' @return A \code{metalog} object with elements
#' \item{params}{A list of the parameters used to create the metalog object}
#' \item{dataValues}{a dataframe with the first column the raw data, second
#'                   column the cumulative probabilities and the third the z
#'                   vector}
#' \item{Y}{The Y matrix values for each quantile and term}
#' \item{A}{a dataframe of coefficients for each metalog distribution}
#' \item{M}{a dataframe of quantiles (M) and probabilities (m) indexed for each
#'          term (i.e. M3,m3 for the third term)}
#' \item{GridPlotCDF()}{a function that displays a grid plot of the CDF for each
#'                      term}
#' \item{VGridPlotPDF()}{a function that displays a gird plot of the PDF for
#'                       each term}
#' \item{Validation}{a vector of yes/no indicators of the valid distributions
#'                   for each term}
#'
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 13)
metalog <- function(x,
                    bounds = c(0,1),
                    boundedness = 'u',
                    term_limit = 13,
                    term_lower_bound = 2,
                    step_len = 0.01,
                    probs = NA,
                    fit_method = 'any') {
  # Input validation
  if (class(x) != 'numeric') {
    stop('Input x must be a numeric vector!')
  }

  if (class(bounds) != 'numeric') {
    stop('Error: bounds must be a numeric vector!')
  }

  if (term_limit %% 1 != 0) {
    stop('Error: term_limit parameter should be an integer between 3 and 30')
  }

  if (term_lower_bound %% 1 != 0) {
    stop('Error: term_lower_bound parameter should be an integer')
  }

  if (length(which(is.na(probs)))==0 &
      (length(probs) != length(x))) {
    stop('Error: probs vector and x vector must be the same length')
  }

  if (length(which(is.na(probs)))==0) {
    if (class(probs) != 'numeric') {
      stop('Error: input probabilites must be a numeric vector!')
    }

    if (max(probs) > 1 | min(probs) < 0) {
      stop('Error: input probabilites must have values between, ',
           'not including, 0 and 1')
    }
  }

  if (length(x) <= 2) {
    stop('Error: input x must be of length 3 or greater')
  }

  if (length(bounds) != 2 & boundedness == 'b') {
    stop('Error: must supply only upper and lower bounds as a numeric vector ',
         '(i.e. c(0,30))')
  }

  if (max(bounds) < min(bounds) & boundedness == 'b') {
    stop('Error: upper bound must be greater than lower bound')
  }

  if (min(x) < min(bounds) & boundedness == 'b') {
    stop('Error: lower bound must be less that the smallest value of x')
  }

  if (max(bounds) < max(x) & boundedness == 'b') {
    stop('Error: upper bound must be greater than the largest value of x')
  }
  if (length(bounds) != 1 &
      (boundedness == 'su' | boundedness == 'sl')) {
    stop('Error: must supply one bound')
  }

  if (boundedness == 'su') {
    bounds <- c(min(x), bounds)
  }

  if (boundedness == 'sl') {
    bounds <- c(bounds, max(x))
  }

  if (boundedness != 'u' &
      boundedness != 'su' &
      boundedness != 'sl' &
      boundedness != 'b') {
    stop('Error: boundedness parameter must be u, su, sl or b only')
  }

  if (max(x) > bounds[2] &
      boundedness == 'su') {
    stop('Error: for semi-upper bounded the upper bound must be greater than ',
         'the largest value in x')
  }

  if (min(x) < bounds[1] &
      boundedness == 'sl') {
    stop('Error: for semi-lower bounded the lower bound must be less than the ',
         'smallest value in x')
  }

  if (term_limit < 3) {
    stop('Error: term_limit should be 3 or greater')
  }

  if (term_limit > 30) {
    stop('Error: term_limit parameter should be less than 30')
  }

  if (term_limit > length(x)) {
    stop('Error: term_limit must be less than or equal to the length of the ',
         'vector x')
  }

  if (term_lower_bound > term_limit) {
    stop('Error: term_lower_bound must be less than or equal to term_limit')
  }

  if (term_lower_bound < 2) {
    stop('Error: term_lower_bound must have a value of 2 or greater')
  }

  if (step_len < 0.001 |
      step_len > 0.01) {
    stop('Error: step_len must be >= to 0.001 and <= to 0.01')
  }

  if (fit_method !='OLS' &
      fit_method !='LP' &
      fit_method !='any') {
    stop('Error: fit_method can only be values OLS, LP or any')
  }

  # Create a list to hold all the objects
  myList <- list()
  myList$params$bounds <- bounds
  myList$params$boundedness <- boundedness
  myList$params$term_limit <- term_limit
  myList$params$term_lower_bound <- term_lower_bound
  myList$params$step_len <- step_len
  myList$params$fit_method <- fit_method
  # Handle the probabilites --- this also converts x as a data frame
  if (length(which(is.na(probs)))!=0) {
    x <- MLprobs(x, step_len = step_len)
  } else{
    x <- as.data.frame(x)
    x$probs <- probs
  }

  # Build the z vector based on the boundedness
  if (boundedness == 'u') {
    x$z <- x[, 1]
  }

  if (boundedness == 'sl') {
    x$z <- log(x[, 1] - bounds[1])
  }

  if (boundedness == 'su') {
    x$z <- (-log(bounds[2] - x[, 1]))
  }

  if (boundedness == 'b') {
    x$z <- log((x[, 1] - bounds[1]) / (bounds[2] - x[, 1]))
  }

  myList$dataValues <- x
  Y <- data.frame(y1 = rep(1, nrow(x)))

  # Construct the Y Matrix initial values
  Y$y2 <- (log(x$probs / (1 - x$probs)))
  Y$y3 <- (x$probs - 0.5) * Y$y2

  if (term_limit > 3) {
    Y$y4 <- x$probs - 0.5
  }

  # Complete the values through the term limit
  if (term_limit > 4) {
    for (i in 5:(term_limit)) {
      y <- paste0('y', i)
      if (i %% 2 != 0) {
        Y[`y`] <- Y$y4 ^ (i %/% 2)
      }
      if (i %% 2 == 0) {
        z <- paste0('y', (i - 1))
        Y[`y`] <- Y$y2 * Y[`z`]
      }
    }
  }
  myList$Y <- Y

  # Build a vectors for each term and
  # build the metalog m(pdf) and M(quantile) dataframes
  myList <- a_vector_OLS_and_LP(
    myList,
    term_limit = term_limit,
    term_lower_bound = term_lower_bound,
    bounds = bounds,
    boundedness = boundedness,
    fit_method = fit_method,
    diff_error = .001,
    diff_step = 0.001
  )

  class(myList) <- append("metalog", class(myList))

  return(myList)
}
