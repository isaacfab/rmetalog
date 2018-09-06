#' Create random samples from a metalog object
#'
#' @param x metalog object created from \code{r_metalog()}
#' @param n number of observations (default is 1)
#' @param term which metalog distribution to sample from
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
#' rmetalog_summary(myMetalog)
rmetalog_sample <- function(x,n=1,term=3) {
  UseMethod("rmetalog_sample",x)
}

#' @export
rmetalog_sample.default <- function(x,n=1,term=3){
  print('Object must be of calss rmetalog')
}

#' @export
rmetalog_sample.rmetalog <- function(x,n=1,term=3){
  print('Object is of calss rmetalog')
  print(n)
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
#' @param x metalog object created from \code{r_metalog()}
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
