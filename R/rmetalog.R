#base rMetalog function here!

#' Fit the Metalog distribution
#'
#' @param x vector of numeric data
#' @param step_len size of steps to evaluate the distribution (between 0 and 1)
#' @param probs (Optional) probability quantiles, same length as \code{x}
#' @param term_limit integer between 5 and 30, specifying number of metalog
#'   terms to evaluate (default: 16)
#' @param bounds numeric vector of size two, indicating lower/upper bounds
#' @param boundedness character string specifying unbounded, semi-bounded upper,
#'   semi-bounded lower or bounded; accepts values \code{u}, \code{su},
#'   \code{sl} and \code{b} (default: 'u')
#'
#' @return A list object with elements
#' \item{Y}{data frame}
#' \item{A}{data frame}
#' \item{M}{data frame}
#' \item{Validation}{character vector (?)}
#'
#' @export
#'
#' @examples
#' # Load example data
#' data("fishSize")
#'
#' # Create a bounded metalog object
#' myMetalog <- rMetalog(fishSize$FishSize,
#'                       step_len = .001,
#'                       bounds=c(0, 60),
#'                       boundedness = 'b',
#'                       term_limit = 16)
rMetalog <-
  function(x,
           probs = 0,
           step_len = .01,
           term_limit = 3,
           bounds = c(0,1),
           boundedness = 'u') {

#create a list to hold all the objects
myList<-list()

################# inital error checking ################
if(class(x)!='numeric'){
  return(print('Error: input x must be a numeric vector!'))
}

if(class(probs)!='numeric'){
  return(print('Error: input probabilites must be a numeric vector!'))
}

if(boundedness!='u'&class(bounds)!='numeric'){
  return(print('Error: bounds must be a numeric vector!'))
}

if(probs[1]!=0&(length(probs)!=length(x))){
  return(print('Error: probability vector and vector x must be the same length'))
}

if(max(probs)>1|min(probs)<0){
  return(print('Error: input probabilites have values between 0 and 1'))
}

if(length(x)<=2){
  return(print('Error: input x must be of length 3 or greater'))
}

if(length(bounds)!=2&boundedness=='b'){
  return(print('Error: must supply upper and lower bounds as a numeric vector (i.e. c(0,30))'))
}

if(length(bounds)!=1&(boundedness=='su'|boundedness=='sl')){
  return(print('Error: must supply one bound'))
}

if(boundedness=='su'){
  bounds<-c(min(x),bounds)
}

if(boundedness=='sl'){
  bounds<-c(bounds,max(x))
}

if(boundedness!='u'&boundedness!='su'&boundedness!='sl'&boundedness!='b'){
  return(print('Error: boundedness parameter must be u, su, sl or b only'))
}

if(max(x)>bounds[2]&boundedness=='su'){
  return(print('Error: for semi-upper bounded the upper bound must be greater than the largest value in x'))
}

if(min(x)<bounds[1]&boundedness=='sl'){
  return(print('Error: for semi-lower bounded the lower bound must be less than the smallest value in x'))
}

if(term_limit%%1!=0){
  return(print('Error: term_limit parameter should be an ineteger between 3 and 30'))
}
if(term_limit<3){
  return(print('Error: term_limit should be 3 or greater'))
}
if(term_limit>30){
  return(print('Error: term_limit parameter should be less than 30'))
}
if(term_limit>length(x)){
  return(print('Error: term_limit must be less than or equal to the length of the vector x'))
}

###############handle the probabilites###############
#this also converts x as a data frame
   if(probs[1] == 0){
     x<-MLprobs(x,step_len=step_len)
   } else{
     x<-as.data.frame(x)
     x$probs<-probs
   }


#########build the z vector based on the boundedness###########
   if(boundedness=='u'){
     x$z<-x[,1]
   }
   if(boundedness=='sl'){
     x$z<-log(x[,1]-bounds[1])
   }
  if(boundedness=='su'){
    x$z<-(-log(bounds[2]-x[,1]))
  }
  if(boundedness=='b'){
    x$z<-log((x[,1]-bounds[1])/(bounds[2]-x[,1]))
  }

################construct the Y Matrix initial values################
pbY<-progress::progress_bar$new(total=(term_limit-1))
  x$y1<-1
  x$y2<-(log(x$probs/(1-x$probs)))
  x$y3<-(x$probs-0.5)*x$y2
  pbY$tick()
if(term_limit>3){
  x$y4<-x$probs-0.5
  pbY$tick()
}
#####complete the values through the term limit#####
if(term_limit>4){
    for (i in 5:(term_limit)){
        pbY$tick()
        y<-paste0('y',i)
        if(i %% 2 != 0){
         x[`y`]<-x$y4^(i%/%2)
        }
        if(i %% 2 == 0){
         z<-paste0('y',(i-1))
         x[`y`]<-x$y2*x[`z`]
        }
    }
}
  myList$Y<-x

###########build a vectors for each term###########
A<-aVectorsMetalogLP(myList$Y,term_limit=term_limit,diff_error=.001,diff_step=0.001)
  myList$A<-A

##############build the metalog m(pdf) and M(quantile) dataframes###############
y<-seq(step_len,(1-step_len),step_len)

Mh<-data.frame()
print('Building distribution functions and samples', row.names=FALSE)
pbP<-progress::progress_bar$new(total=(term_limit-1))
for(i in 2:term_limit){
  pbP$tick()
  a_name<-paste0('a',i)
  m_name<-paste0('m',i)
  M_name<-paste0('M',i)

#build pdf
  m<-pdfMetalog(myList$A[`a_name`][,1],y[1],term_limit,bounds=bounds,boundedness=boundedness)

  for(j in 2:length(y)){
    temp<-pdfMetalog(myList$A[`a_name`][,1],y[j],term_limit,bounds=bounds,boundedness=boundedness)
    m<-c(m,temp)
  }

#build quantile values
  M<-quantileMetalog(myList$A[`a_name`][,1],y[1],term_limit,bounds=bounds,boundedness=boundedness)

  for(j in 2:length(y)){
    temp<-quantileMetalog(myList$A[`a_name`][,1],y[j],term_limit,bounds=bounds,boundedness=boundedness)
    M<-c(M,temp)
  }


#add traling and leading zero's for pdf bounds
  if(boundedness=='sl'){
    m<-c(0,m)
    M<-c(bounds[1],M)
  }
  if(boundedness=='su'){
    m<-c(m,0)
    M<-c(M,bounds[2])
  }
  if(boundedness=='b'){
    m<-c(0,m,0)
    M<-c(bounds[1],M,bounds[2])
  }

  if(length(Mh)==0){
    Mh<-data.frame(m2=m,M2=M)
  }

  if(length(M)!=0){
     Mh[`m_name`]<-m
     Mh[`M_name`]<-M
  }
}#close for loop

#adding the y values for the bounded models
if(boundedness=='sl'){
  y<-c(0,y)
}
if(boundedness=='su'){
  y<-c(y,1)
}
if(boundedness=='b'){
  y<-c(0,y,1)
}


Mh$y<-y

#return values
myList$M<-Mh

#build plots
InitalResults<-data.frame(term=(rep(c('2 Terms'),length(Mh[,1]))),pdfValues=Mh$m2,quantileValues=Mh$M2,cumValue=Mh$y)

for(i in 2:(length(Mh[1,]-1)/2)){
  TempResults<-data.frame(term=(rep(paste0((i+1),' Terms'),length(Mh[,1]))),pdfValues=Mh[,(i*2-1)],quantileValues=Mh[,(i*2)],cumValue=Mh$y)
  InitalResults<-rbind(InitalResults,TempResults)
}

# The base plot
q <- ggplot2::ggplot(InitalResults, ggplot2::aes(x=quantileValues, y=pdfValues)) + ggplot2::geom_line()
#q$term<-as.factor(q$term)

# Faceted using subpanels
q<-q + ggplot2::facet_wrap(~term,ncol=4,scales="free_y")

q

myList$GridPlotPDF<-q

# The base plot
q <- ggplot2::ggplot(InitalResults, ggplot2::aes(x=quantileValues, y=cumValue)) + ggplot2::geom_line()
#q$term<-as.factor(q$term)

# Faceted using subpanels
q<-q + ggplot2::facet_wrap(~term,ncol=4,scales="free_y")

q

myList$GridPlotCDF<-q
#########pdf validation################
y<-c()
t<-c()
for(i in 2:term_limit){
  m<-paste0('m',i)
  x<-pdfMetalogValidation(myList$M[`m`])
  t<-c(t,i)
  y<-c(y,x)

}
val<-data.frame(term=t,valid=y)
myList$Validation<-val


return(myList)
}


