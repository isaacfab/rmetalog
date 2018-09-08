#' Create random samples from an rmetalog distribution object
#'
#' @param my_metalog metalog object created from \code{r_metalog()}
#' @param n number of observations (default is 1)
#' @param term which metalog distribution to sample from
#'
#' @return A numeric vector of n random samples from a selected distribution
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- r_metalog(fishSize$FishSize,
#'                       bounds=c(0, 60),
#'                       boundedness = 'b',
#'                       term_limit = 9,
#'                       term_lower_bound = 9)
#'
#' s <- r_metalog_sample(myMetalog, n=1000, term = 9)
#' hist(s)
r_metalog_sample <- function(my_metalog,n=1,term=3) {
  UseMethod("r_metalog_sample",my_metalog)
}

#' @export
r_metalog_sample.default <- function(my_metalog,n=1,term=3){
  print('Object must be of calss rmetalog')
}

#' @export
r_metalog_sample.rmetalog <- function(my_metalog,n=1,term=3){
  ###################some error checking######################
  valid_terms<-my_metalog$Validation$term
  if(class(n)!='numeric'|n<1|n%%1!=0){
    stop('Error: n must be a positive numeric interger')
  }
  if(class(term)!='numeric'|term<2|term%%1!=0|!(term %in% valid_terms)|length(term)>1){
    stop(paste0('Error: term must be a single positive numeric interger contained in the metalog object. Available terms are: ',
                valid_terms))
  }
  x<-runif(n)
  Y<-data.frame(y1=rep(1,n))
  ################construct the Y Matrix initial values################

  Y$y2<-(log(x/(1-x)))
  Y$y3<-(x-0.5)*Y$y2

  if(term>3){
    Y$y4<-x-0.5
  }
  #####complete the values through the term limit#####
  if(term>4){
    for (i in 5:(term)){

      y<-paste0('y',i)
      if(i %% 2 != 0){
        Y[`y`]<-Y$y4^(i%/%2)
      }
      if(i %% 2 == 0){
        z<-paste0('y',(i-1))
        Y[`y`]<-Y$y2*Y[`z`]
      }
    }
  }
 Y<-as.matrix(Y)
 amat<-paste0('a',term)
 a<-as.matrix(my_metalog$A[`amat`])
 s<-Y %*% a
 if(my_metalog$params$boundedness=='sl'){
      s<-my_metalog$params$bounds[1] + exp(s)
 }
 if(my_metalog$params$boundedness=='su'){
      s<-my_metalog$params$bounds[2] - exp(-(s))
 }
 if(my_metalog$params$boundedness=='b'){
     s <- (my_metalog$params$bounds[1]+(my_metalog$params$bounds[2])*exp(s))/(1+exp(s))
 }

 return(s)
}

#' Generate quantiles from a probability from a metalog object
#'
#' @param my_metalog metalog object created from \code{r_metalog()}
#' @param y  y vector of probabilities
#' @param term which metalog distribution to sample from
#'
#' @return A numeric vector of quantiles corresponding to the y probability vector
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- r_metalog(fishSize$FishSize,
#'                       bounds=c(0, 60),
#'                       boundedness = 'b',
#'                       term_limit = 9,
#'                       term_lower_bound = 9)
#'
#' s<-r_metalog_quantile(myMetalog,y=c(0.25,0.5,0.7),term = 9)
r_metalog_quantile <- function(my_metalog,y,term=3) {
  UseMethod("r_metalog_quantile",my_metalog)
}

#' @export
r_metalog_quantile.default <- function(my_metalog,y,term=3){
  print('Object must be of class rmetalog')
}

#' @export
r_metalog_quantile.rmetalog <- function(my_metalog,y,term=3){
  ###################some error checking######################
  valid_terms<-my_metalog$Validation$term
  if(class(y)!='numeric'|max(y)>=1|min(y)<=0){
    stop('Error: y must be a positive numeric vector between 0 and 1')
  }
  if(class(term)!='numeric'|term<2|term%%1!=0|!(term %in% valid_terms)|length(term)>1){
    stop(paste0('Error: term must be a single positive numeric interger contained in the metalog object. Available terms are: ',
                valid_terms))
  }

  Y<-data.frame(y1=rep(1,length(y)))
  ################construct the Y Matrix initial values################

  Y$y2<-(log(y/(1-y)))
  Y$y3<-(y-0.5)*Y$y2

  if(term>3){
    Y$y4<-(y-0.5)
  }
  #####complete the values through the term limit#####
  if(term>4){
    for (i in 5:(term)){

      y<-paste0('y',i)
      if(i %% 2 != 0){
        Y[`y`]<-Y$y4^(i%/%2)
      }
      if(i %% 2 == 0){
        z<-paste0('y',(i-1))
        Y[`y`]<-Y$y2*Y[`z`]
      }
    }
  }
  Y<-as.matrix(Y)
  amat<-paste0('a',term)
  a<-as.matrix(my_metalog$A[`amat`])
  s<-Y %*% a
  if(my_metalog$params$boundedness=='sl'){
    s<-my_metalog$params$bounds[1] + exp(s)
  }
  if(my_metalog$params$boundedness=='su'){
    s<-my_metalog$params$bounds[2] - exp(-(s))
  }
  if(my_metalog$params$boundedness=='b'){
    s <- (my_metalog$params$bounds[1]+(my_metalog$params$bounds[2])*exp(s))/(1+exp(s))
  }

  return(s)
}


#' Summary of the metalog object
#'
#' @param object metalog object created from \code{r_metalog()}
#' @param ... other stuff
#'
#' @return A summary of the object
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- r_metalog(fishSize$FishSize,
#'                        bounds=c(0, 60),
#'                        boundedness = 'b',
#'                        term_limit = 13)
#'
#' summary(myMetalog)
summary.rmetalog <- function(object, ...) {
  cat(' -----------------------------------------------\n',
      'Summary of Metalog Distribution Object\n',
      '-----------------------------------------------\n',
      '\nParameters\n',
      'Term Limit: ', object$params$term_limit, '\n',
      'Term Lower Bound: ', object$params$term_lower_bound, '\n',
      'Boundedness: ', object$params$boundedness, '\n',
      'Bounds (only used based on boundedness): ', object$params$bounds, '\n',
      'Step Length for Distribution Summary: ', object$params$step_len, '\n',
      '\n\n Validation and Fit Method\n'
  )
  print(object$Validation, row.names = FALSE)
}

#' Plot of the metalog object
#'
#' @param x metalog object created using \code{r_metalog()}
#' @param ... ignored; included for S3 generic/method consistency
#'
#' @return A summary plot of the CDF and PDF for each term
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- r_metalog(fishSize$FishSize,
#'                        bounds=c(0, 60),
#'                        boundedness = 'b',
#'                        term_limit = 13)
#'
#' r_metalog_plot(myMetalog)
plot.rmetalog <- function(x, ...){
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
    for (i in 2:(length(x$M[1, ] - 1) / 2)) {
      TempResults <-
        data.frame(
          term = (rep(paste0((x$params$term_lower_bound + (i - 1)), ' Terms'
          ), length(x$M[, 1]))),
          pdfValues = x$M[, (i * 2 - 1)],
          quantileValues = x$M[, (i * 2)],
          cumValue = x$M$y
        )

      InitalResults <- rbind(InitalResults, TempResults)
    }
  }

  # The base plot
  p <-
    ggplot2::ggplot(InitalResults, ggplot2::aes(x = quantileValues, y = pdfValues)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::xlab("Quantile Values") +
    ggplot2::ylab("PDF Values") +
    ggplot2::theme_bw()

  # Faceted using subpanels
  p <- p + ggplot2::facet_wrap( ~ term, ncol = 4, scales = "free_y")

  # readline(prompt="Press [enter] to see CDF plot")

  # The base plot
  q <-
    ggplot2::ggplot(InitalResults, ggplot2::aes(x = quantileValues, y = cumValue)) +
    ggplot2::geom_line(colour = "blue") +
    ggplot2::xlab("Quantile Values") +
    ggplot2::ylab("CDF Values") +
    ggplot2::theme_bw()
  #q$term<-as.factor(q$term)

  # Faceted using subpanels
  q <- q + ggplot2::facet_wrap( ~ term, ncol = 4, scales = "free_y")

  list(pdf = p, cdf = q)
}
