# Given an a vector, produce a dataframe with pdf and quantile values for cumulants
pdf_quantile_builder <-
  function(temp, y, term_limit, bounds, boundedness) {
    myList <- list()

    # Build pdf
    m <- pdfMetalog(temp, y[1], term_limit, bounds = bounds,
                    boundedness = boundedness)

    for (j in 2:length(y)) {
      tempPDF <- pdfMetalog(temp, y[j], term_limit, bounds = bounds,
                            boundedness = boundedness)
      m <- c(m, tempPDF)
    }

    # Build quantile values
    M <- quantileMetalog(temp, y[1], term_limit, bounds = bounds,
                         boundedness = boundedness)

    for (j in 2:length(y)) {
      tempQant <- quantileMetalog(temp, y[j], term_limit, bounds = bounds,
                                  boundedness = boundedness)
      M <- c(M, tempQant)
    }

    # Add trailing and leading zero's for pdf bounds
    if (boundedness == 'sl') {
      m <- c(0, m)
      M <- c(bounds[1], M)
    }

    if (boundedness == 'su') {
      m <- c(m, 0)
      M <- c(M, bounds[2])
    }

    if (boundedness == 'b') {
      m <- c(0, m, 0)
      M <- c(bounds[1], M, bounds[2])
    }

    # Add y values for bounded models
    if (boundedness == 'sl') {
      y <- c(0, y)
    }

    if (boundedness == 'su') {
      y <- c(y, 1)
    }

    if (boundedness == 'b') {
      y <- c(0, y, 1)
    }

    myList$m <- m
    myList$M <- M
    myList$y <- y

    # PDF validation
    myList$valid <- pdfMetalogValidation(myList$m)

    return(myList)
  }

# pdf validation function
# call this feasibility
pdfMetalogValidation <- function(x) {
  y <- min(x)
  if (y >= 0) {
    return('yes')
  }
  if (y < 0) {
    return('no')
  }
}
