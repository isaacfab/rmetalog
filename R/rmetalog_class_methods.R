#' Create random samples from an rmetalog distribution object
#'
#' @param x metalog object created from \code{r_metalog()}
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
#'                       term_limit = 13)
#'
#' s<-rmetalog_sample(myMetalog,1000,9)
#' hist(s)
rmetalog_sample <- function(x,n=1,term=3) {
  UseMethod("rmetalog_sample",x)
}

#' @export
rmetalog_sample.default <- function(x,n=1,term=3){
  print('Object must be of calss rmetalog')
}

#' @export
rmetalog_sample.rmetalog <- function(myList,n=1,term=3){
  ###################some error checking######################
  valid_terms<-myList$Validation$term
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
 a<-as.matrix(myList$A[`amat`])
 s<-Y %*% a
 if(myList$params$boundedness=='sl'){
      s<-my_metalog$params$bounds[1] + exp(s)
 }
 if(myList$params$boundedness=='su'){
      s<-my_metalog$params$bounds[2] - exp(-(s))
 }
 if(myList$params$boundedness=='b'){
     s <- (my_metalog$params$bounds[1]+(my_metalog$params$bounds[2])*exp(s))/(1+exp(s))
 }

 return(s)
}

#' Generate quantiles from a probability from a metalog object
#'
#' @param x metalog object created from \code{r_metalog()}
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
#'                       term_limit = 13)
#'
#' s<-rmetalog_quantile(myMetalog)
rmetalog_quantile <- function(x,y,term=3) {
  UseMethod("rmetalog_quantile",x)
}

#' @export
rmetalog_quantile.default <- function(x,y,term=3){
  print('Object must be of calss rmetalog')
}

#' @export
rmetalog_quantile.rmetalog <- function(myList,y,term=3){
  ###################some error checking######################
  valid_terms<-myList$Validation$term
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
  a<-as.matrix(myList$A[`amat`])
  s<-Y %*% a
  if(myList$params$boundedness=='sl'){
    s<-my_metalog$params$bounds[1] + exp(s)
  }
  if(myList$params$boundedness=='su'){
    s<-my_metalog$params$bounds[2] - exp(-(s))
  }
  if(myList$params$boundedness=='b'){
    s <- (my_metalog$params$bounds[1]+(my_metalog$params$bounds[2])*exp(s))/(1+exp(s))
  }

  return(s)
}


#' Summary of the metalog object
#'
#' @param x metalog object created from \code{r_metalog()}
#'
#' @return A summary of the object and what it contains
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
#'                       term_limit = 13)
#'
#' summary(myMetalog)
summary.rmetalog <- function(x){
  cat(' -----------------------------------------------\n',
      'Summary of Metalog Distribution Object\n',
      '-----------------------------------------------\n',
      '\nParameters\n',
      'Term Limit: ',x$params$term_limit, '\n',
      'Term Lower Bound: ',x$params$term_lower_bound, '\n',
      'Boundedness: ',x$params$boundedness, '\n',
      'Bounds (only used based on boundedness): ',x$params$bounds, '\n',
      'Step Length for Distribution Summary: ',x$params$step_len, '\n',
      '\n\n Validation and Fit Method\n'
  )
  print(x$Validation,row.names=FALSE)

}

#' Plot of the metalog object
#'
#' @param x metalog object created using \code{r_metalog()}
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
#'                       bounds=c(0, 60),
#'                       boundedness = 'b',
#'                       term_limit = 13)
#'
#' plot(myMetalog)
plot.rmetalog <- function(myList){
  #build plots
  InitalResults<-data.frame(term=(rep(paste0(myList$params$term_lower_bound,' Terms'),length(myList$M[,1]))),
                            pdfValues=myList$M[,1],
                            quantileValues=myList$M[,2],
                            cumValue=myList$M$y)
  if(ncol(myList$M)>3){
    for(i in 2:(length(myList$M[1,]-1)/2)){
      TempResults<-data.frame(term=(rep(paste0((myList$params$term_lower_bound+(i-1)),' Terms'),length(myList$M[,1]))),
                              pdfValues=myList$M[,(i*2-1)],
                              quantileValues=myList$M[,(i*2)],
                              cumValue=myList$M$y)

      InitalResults<-rbind(InitalResults,TempResults)
    }
  }
  # The base plot
  q <- ggplot2::ggplot(InitalResults, ggplot2::aes(x=quantileValues, y=pdfValues)) +
    ggplot2::geom_line(colour="blue") +
    ggplot2::xlab("Quantile Values") +
    ggplot2::ylab("PDF Values") +
    ggplot2::theme_bw()

  # Faceted using subpanels
  q <- q + ggplot2::facet_wrap(~term,ncol=4,scales="free_y")

  print(q)

  readline(prompt="Press [enter] to see CDF plot")

  # The base plot
  q <- ggplot2::ggplot(InitalResults, ggplot2::aes(x=quantileValues, y=cumValue)) +
    ggplot2::geom_line(colour="blue") +
    ggplot2::xlab("Quantile Values") +
    ggplot2::ylab("CDF Values") +
    ggplot2::theme_bw()
  #q$term<-as.factor(q$term)

  # Faceted using subpanels
  q <- q + ggplot2::facet_wrap(~term,ncol=4,scales="free_y")

  print(q)

}
