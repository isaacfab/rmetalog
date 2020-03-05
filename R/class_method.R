#' Create random samples from an metalog distribution object
#'
#' @param m metalog object created from \code{metalog()}
#' @param n number of observations (default is 1)
#' @param term which metalog distribution to sample from
#'
#' @return A numeric vector of n random samples from a selected distribution
#'
#' @export
#'
#' @examples
#' # Load example data
#' \dontrun{
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 9,
#'                      term_lower_bound = 9)
#'
#' s <- rmetalog(myMetalog, n=1000, term = 9)
#' hist(s)
#' }
rmetalog <- function(m, n = 1, term = 3) {
  UseMethod("rmetalog", m)
}

#' @export
rmetalog.default <- function(m, n = 1, term = 3) {
  print('Object must be of calss metalog')
}

#' @export
rmetalog.metalog <- function(m, n = 1, term = 3){
  # Input validation
  valid_terms <- m$Validation$term
  valid_terms_printout <- paste(valid_terms, collapse = " ")
  if (class(n) != 'numeric' | n < 1 | n %% 1 != 0) {
    stop('Error: n must be a positive numeric interger')
  }
  if (class(term) != 'numeric' |
      term < 2 | term %% 1 != 0 | !(term %in% valid_terms) |
      length(term) > 1) {
    stop(
      paste('Error: term must be a single positive numeric interger contained',
            'in the metalog object. Available terms are:', valid_terms_printout)
    )
  }
  x <- stats::runif(n)
  Y <- data.frame(y1 = rep(1, n))

  # Construct initial Y Matrix values
  Y$y2 <- (log(x / (1 - x)))

  if (term > 2) {
    Y$y3 <- (x - 0.5) * Y$y2
  }

  if (term > 3) {
    Y$y4 <- x - 0.5
  }

  # Complete the values through the term limit
  if (term > 4) {
    for (i in 5:(term)) {
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

  Y <- as.matrix(Y)
  amat <- paste0('a', term)
  a <- as.matrix(m$A[`amat`])
  s <- Y %*% a[1:term]

  if (m$params$boundedness == 'sl') {
    s <- m$params$bounds[1] + exp(s)
  }

  if (m$params$boundedness == 'su') {
    s <- m$params$bounds[2] - exp(-(s))
  }

  if (m$params$boundedness == 'b') {
    s <-
      (m$params$bounds[1] + (m$params$bounds[2]) * exp(s)) /
      (1 + exp(s))
  }

  return(as.numeric(s))
}


#' Generate quantiles with a probability from a metalog object
#'
#' @param m metalog object created from \code{metalog()}
#' @param y  vector of probabilities
#' @param term which metalog distribution to sample from
#'
#' @return A numeric vector of quantiles corresponding to the y probability
#'   vector
#'
#' @export
#'
#' @examples
#' # Load example data
#' \dontrun{
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 9,
#'                      term_lower_bound = 9)
#'
#' s <- qmetalog(myMetalog,y=c(0.25,0.5,0.7),term = 9)
#' }
qmetalog <- function(m, y, term = 3) {
  UseMethod("qmetalog", m)
}

#' @export
qmetalog.default <- function(m, y, term = 3){
  print('Object must be of class metalog')
}

#' @export
qmetalog.metalog <- function(m, y, term = 3){
  # Input validation
  valid_terms <- m$Validation$term
  valid_terms_printout <- paste(valid_terms, collapse = " ")
  if (class(y) != 'numeric' | max(y) >= 1 | min(y) <= 0) {
    stop('Error: y must be a positive numeric vector between 0 and 1')
  }

  if (class(term) != 'numeric' |
      term < 2 | term %% 1 != 0 | !(term %in% valid_terms) |
      length(term) > 1) {
    stop(
      paste('Error: term must be a single positive numeric interger contained',
            'in the metalog object. Available terms are:', valid_terms_printout)
    )
  }

  Y <- data.frame(y1 = rep(1, length(y)))

  # Construct the Y Matrix initial values
  Y$y2 <- (log(y / (1 - y)))

  if (term > 2) {
    Y$y3 <- (y - 0.5) * Y$y2
  }

  if (term > 3) {
    Y$y4 <- (y - 0.5)
  }

  # Complete the values through the term limit
  if (term > 4) {
    for (i in 5:(term)) {
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

  Y <- as.matrix(Y)
  amat <- paste0('a', term)
  a <- as.matrix(m$A[`amat`])
  s <- Y %*% a[1:term]

  if (m$params$boundedness == 'sl') {
    s <- m$params$bounds[1] + exp(s)
  }

  if (m$params$boundedness == 'su') {
    s <- m$params$bounds[2] - exp(-(s))
  }

  if (m$params$boundedness == 'b') {
    s <- (m$params$bounds[1] + (m$params$bounds[2]) * exp(s)) / (1 + exp(s))
  }

  return(as.numeric(s))
}

#' Generate probabilities with quantiles from a metalog object.
#' This is done through a newtons method approximation.
#'
#' @param m metalog object created from \code{metalog()}
#' @param q  vector of quantiles
#' @param term which metalog distribution to sample from
#'
#' @return A numeric vector of probabilities corresponding to the q quantile
#'   vector
#'
#' @export
#'
#' @examples
#' # Load example data
#' \dontrun{
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 9,
#'                      term_lower_bound = 9)
#'
#' s <- pmetalog(myMetalog,q=c(3,10,25),term = 9)
#' }
pmetalog <- function(m, q, term = 3) {
  UseMethod("pmetalog", m)
}

pmetalog.default <- function(m, q, term = 3){
  print('Object must be of class metalog')
}

#' @export
pmetalog.metalog <- function(m, q, term = 3){
  # Input validation
  valid_terms <- m$Validation$term
  if (class(q) != 'numeric') {
    stop('Error: q must be a positive numeric vector between 0 and 1')
  }

  if (class(term) != 'numeric' |
      term < 2 | term %% 1 != 0 | !(term %in% valid_terms) |
      length(term) > 1) {
    stop(
      cat('Error: term must be a single positive numeric interger contained',
            'in the metalog object. Available terms are:',
            valid_terms)
    )
  }

 qs<-sapply(q,newtons_method_metalog,m=m,t=term)
 return(qs)
}

#' Generate density values with quantiles from a metalog object.
#' This is done through a newtons method approximation.
#'
#' @param m metalog object created from \code{metalog()}
#' @param q  y vector of quantiles
#' @param term which metalog distribution to sample from
#'
#' @return A numeric vector of probabilities corresponding to the q quantile
#'   vector
#'
#' @export
#'
#' @examples
#' # Load example data
#' \dontrun{
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 9,
#'                      term_lower_bound = 9)
#'
#' s <- dmetalog(myMetalog,q=c(3,10,25),term = 9)
#' }
dmetalog <- function(m, q, term = 3) {
  UseMethod("dmetalog", m)
}

dmetalog.default <- function(m, q, term = 3){
  print('Object must be of class metalog')
}

#' @export
dmetalog.metalog <- function(m, q, term = 3){
  # Input validation
  valid_terms <- m$Validation$term
  if (class(q) != 'numeric') {
    stop('Error: q must be a numeric vector')
  }

  if (class(term) != 'numeric' |
      term < 2 | term %% 1 != 0 | !(term %in% valid_terms) |
      length(term) > 1) {
    stop(
      paste('Error: term must be a single positive numeric interger contained',
            'in the metalog object. Available terms are:',
            valid_terms)
    )
  }

  qs<-sapply(q,newtons_method_metalog,m=m,term=term)
  ds<-sapply(qs,pdfMetalog_density, m=m,t=term)
  return(ds)
}

#' Summary of the metalog object
#'
#' @param object metalog object created from \code{metalog()}
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return A summary of the object
#'
#' @export
#'
#' @examples
#' # Load example data
#' \dontrun{
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 13)
#'
#' summary(myMetalog)
#' }
summary.metalog <- function(object, ...) {
  cat(' -----------------------------------------------\n',
      'Summary of Metalog Distribution Object\n',
      '-----------------------------------------------\n',
      '\nParameters\n',
      'Term Limit: ', object$params$term_limit, '\n',
      'Term Lower Bound: ', object$params$term_lower_bound, '\n',
      'Boundedness: ', object$params$boundedness, '\n',
      'Bounds (only used based on boundedness): ', object$params$bounds, '\n',
      'Step Length for Distribution Summary: ', object$params$step_len, '\n',
      'Method Use for Fitting: ', object$params$fit_method, '\n',
      'Number of Data Points Used: ', object$params$number_of_data, '\n',
      'Original Data Saved: ', object$params$save_data, '\n',
      '\n\n Validation and Fit Method\n'
  )
  print(object$Validation, row.names = FALSE)
}


#' Plot of the metalog object
#'
#' @param x metalog object created using \code{metalog()}
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return A summary plot of the CDF and PDF for each term
#'
#' @export
#'
#' @examples
#' # Load example data
#' \dontrun{
#' data("fishSize")
#'
#' # Create a bounded metalog object
#'
#' myMetalog <- metalog(fishSize$FishSize,
#'                      bounds=c(0, 60),
#'                      boundedness = 'b',
#'                      term_limit = 13)
#'
#' plot(myMetalog)
#' }
plot.metalog <- function(x, ...) {
  #build plots
  InitalResults <-
    data.frame(
      term = (rep(
        paste0(x$params$term_lower_bound, ' Terms'), length(x$M[, 1])
      )),
      pdfValues = x$M[, 1],
      quantileValues = x$M[, 2],
      cumValue = x$M$y
    )

  if (ncol(x$M) > 3) {
    for (i in 2:(length(x$M[1,] - 1) / 2)) {
      TempResults <-
        data.frame(
          term = (rep(paste((x$params$term_lower_bound + (i - 1)), 'Terms'),
                      length(x$M[, 1]))),
          pdfValues = x$M[, (i * 2 - 1)],
          quantileValues = x$M[, (i * 2)],
          cumValue = x$M$y
        )

      InitalResults <- rbind(InitalResults, TempResults)
    }
  }

  # PDF plot
  p <-
    ggplot2::ggplot(InitalResults, aes(x = quantileValues, y = pdfValues)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::xlab("Quantile Values") +
    ggplot2::ylab("PDF Values") +
    ggplot2::facet_wrap(~ term, ncol = 4, scales = "free_y")

  # The base plot
  q <-
    ggplot2::ggplot(InitalResults, aes(x = quantileValues, y = cumValue)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::xlab("Quantile Values") +
    ggplot2::ylab("CDF Values") +
    ggplot2::facet_wrap(~ term, ncol = 4, scales = "free_y")

  list(pdf = p, cdf = q)
}

