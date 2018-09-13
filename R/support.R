# Supporting functions called inside the metalog function call

# Build the quantiles through a base function
MLprobs <- function(x, step_len) {
  if (class(x) != 'numeric') {
    return(print('Error: input must be a numeric vector!'))
  }
  l <- length(x)
  x <- as.data.frame(x)

  # Need to sort the dataframe for the quantiles to work out
  x <- x[order(x), ]
  x <- as.data.frame(x)

  # Calculate the liklihood as an interpolation
  x$probs <- 0
  for (i in 1:l) {
    if (i == 1) {
      x$probs[i] <- (0.5 / l)
    }
    else{
      x$probs[i] <- (x$probs[i - 1] + (1 / l))
    }
  }

  # If the data is very long we down convert to a smaller but representative
  # vector using the step_len default is 0.01 which is a 109 element vector with
  # fine values in the tail (tailstep)
  if (nrow(x) > 100) {
    y <- seq(step_len, (1 - step_len), step_len)
    tailstep <- (step_len / 10)
    y <- c(seq(tailstep, (min(y) - tailstep), tailstep),
           y,
           seq((max(y) + tailstep), (max(y) + tailstep * 9), tailstep))
    x_new <- stats::quantile(x[, 1], probs = y)
    x <- as.data.frame(x_new)
    x$probs <- y
  }

  return(x)
}

pdfMetalog <- function(a,
                       y,
                       t,
                       bounds = c(),
                       boundedness = 'u') {
  d <- y * (1 - y)
  f <- (y - 0.5)
  l <- log(y / (1 - y))

  # Initiate pdf

  # For the first three terms
  x <- (a[2] / d)
  if (a[3] != 0) {
    x <- x + a[3] * ((f / d) + l)
  }

  # For the fourth term
  if (t > 3) {
    x <- x + a[4]
  }

  # Initalize some counting variables
  e <- 1
  o <- 1

  # For all other terms greater than 4
  if (t > 4) {
    for (i in 5:t) {
      if (i %% 2 != 0) {
        # iff odd
        x <- x + ((o + 1) * a[i] * f^o)
        o <- o + 1
      }

      if (i %% 2 == 0) {
        # iff even
        x <- x + a[i] * (((f^(e + 1)) / d) + (e + 1) * (f^e) * l)
        e <- e + 1
      }
    }
  }

  # Some change of variables here for boundedness

  x <- (x^(-1))

  if (boundedness != 'u') {
    M <- quantileMetalog(a, y, t, bounds = bounds, boundedness = 'u')
  }

  if (boundedness == 'sl') {
    x <- x * exp(-M)
  }

  if (boundedness == 'su') {
    x <- x * exp(M)
  }

  if (boundedness == 'b') {
    x <- (x * (1 + exp(M))^2) / ((bounds[2] - bounds[1]) * exp(M))
  }

  return(x)
}

# Quantile function
quantileMetalog <- function(a,
                            y,
                            t,
                            bounds = c(),
                            boundedness = 'u') {
  # Some values for calculation
  f <- (y - 0.5)
  l <- log(y / (1 - y))

  # For the first three terms
  x <- a[1] + a[2] * l + a[3] * f * l

  # For the fourth term
  if (t > 3) {
    x <- x + a[4] * f
  }

  # Some tracking variables
  o <- 2
  e <- 2

  # For all other terms greater than 4
  if (t > 4) {
    for (i in 5:t) {
      if (i %% 2 == 0) {
        x <- x + a[i] * f^e * l
        e <- e + 1
      }

      if (i %% 2 != 0) {
        x <- x + a[i] * f^o
        o <- o + 1
      }
    }
  }

  if (boundedness == 'sl') {
    x <- bounds[1] + exp(x)
  }

  if (boundedness == 'su') {
    x <- bounds[2] - exp(-x)
  }

  if (boundedness == 'b') {
    x <- (bounds[1] + bounds[2] * exp(x)) / (1 + exp(x))
  }

  return(x)
}

# Function for returning the matrix of differentiation terms
diffMatMetalog <- function(term_limit, step_len) {
  y <- seq(step_len, (1 - step_len), step_len)
  Diff <- data.frame()

  for (i in 1:length(y)) {
    d <- y[i] * (1 - y[i])
    f <- (y[i] - 0.5)
    l <- log(y[i] / (1 - y[i]))

    # Initiate pdf
    diffVector <- 0

    # For the first three terms
    x <- (1 / d)
    diffVector <- c(diffVector, x)

    if (term_limit > 2) {
      diffVector <- c(diffVector, ((f / d) + l))
    }

    # For the fourth term
    if (term_limit > 3) {
      diffVector <- c(diffVector, 1)
    }

    # Initalize some counting variables
    e <- 1
    o <- 1

    # For all other terms greater than 4
    if (term_limit > 4) {
      for (i in 5:term_limit) {
        if (i %% 2 != 0) {
          # iff odd
          diffVector <- c(diffVector, ((o + 1) * f^o))
          o <- o + 1

        }
        if (i %% 2 == 0) {
          # iff even
          diffVector <- c(diffVector,
                          (((f^(e + 1)) / d) + (e + 1) * (f^e) * l))
          e <- e + 1
        }
      }
    }

    Diff <- rbind(Diff, diffVector)
  }

  Diff <- as.matrix(Diff)
  Diff_neg = -(Diff)
  new_Diff <- matrix(c(Diff[, 1], Diff_neg[, 1]), ncol = 2)

  for (c in 2:length(Diff[1, ])) {
    new_Diff <- cbind(new_Diff, Diff[, c])
    new_Diff <- cbind(new_Diff, Diff_neg[, c])
  }

  return(new_Diff)
}

newtons_method_metalog <- function(m,q,term){
  #a simple newtons method application
  alpha_step<-0.01
  err<-0.0000001
  temp_err<-0.1
  y_now<-0.5
  avec<-paste0('a',term)
  a<-m$A[,`avec`]
  i<-1
  while(temp_err>err){
    frist_function<-(quantileMetalog(a,y_now,term,m$params$bounds,m$params$boundedness)-q)
    derv_function<-pdfMetalog(a,y_now,term,m$params$bounds,m$params$boundedness)
    y_next<-y_now-alpha_step*(frist_function*derv_function)
    temp_err<-abs((y_next-y_now))
    if(y_next>1){
      y_next<-0.99999
    }
    if(y_next<0){
      y_next<-0.000001
    }
    y_now<-y_next
    i<-i+1
    if(i>10000){
      stop(paste0('Approximation taking too long, quantile value: ',q,' is to far from distribution median. Try plot() to see distribution.'))
    }
  }

  return(y_now)
}

pdfMetalog_density <- function(m,t,y) {

  avec<-paste0('a',t)
  a<-m$A[,`avec`]
  bounds<-m$params$bounds
  boundedness<-m$params$boundedness

  d <- y * (1 - y)
  f <- (y - 0.5)
  l <- log(y / (1 - y))

  # Initiate pdf

  # For the first three terms
  x <- (a[2] / d)
  if (a[3] != 0) {
    x <- x + a[3] * ((f / d) + l)
  }

  # For the fourth term
  if (t > 3) {
    x <- x + a[4]
  }

  # Initalize some counting variables
  e <- 1
  o <- 1

  # For all other terms greater than 4
  if (t > 4) {
    for (i in 5:t) {
      if (i %% 2 != 0) {
        # iff odd
        x <- x + ((o + 1) * a[i] * f^o)
        o <- o + 1
      }

      if (i %% 2 == 0) {
        # iff even
        x <- x + a[i] * (((f^(e + 1)) / d) + (e + 1) * (f^e) * l)
        e <- e + 1
      }
    }
  }

  # Some change of variables here for boundedness

  x <- (x^(-1))

  if (boundedness != 'u') {
    M <- quantileMetalog(a, y, t, bounds = bounds, boundedness = 'u')
  }

  if (boundedness == 'sl') {
    x <- x * exp(-M)
  }

  if (boundedness == 'su') {
    x <- x * exp(M)
  }

  if (boundedness == 'b') {
    x <- (x * (1 + exp(M))^2) / ((bounds[2] - bounds[1]) * exp(M))
  }

  return(x)
}

