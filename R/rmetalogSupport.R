#Supporting functions called inside the metalog function call

#build the quantiles through a base function
MLprobs <- function(x) {
  if(class(x)!='numeric'){
    return(print('Error: input must be a numeric vector!'))
  }
  l <- length(x)
  x<-as.data.frame(x)
  x$probs <- 0
  for(i in 1:l) {
    if(i==1){
      x$probs[i]<-(0.5/l)
    }
    else{
      x$probs[i]<-(x$probs[i-1]+(1/l))
    }
  }
  return(x)
}

pdfMetalog<-function(a,y,t,bounds=c(),boundedness='u'){
  #error check that a is a numeric vector, y is a number between 0,1 and t is greater than a*
  #some values for calculation
  d<-y*(1-y)
  f<-(y-0.5)
  l<-log(y/(1-y))

  #initiate pdf

  #for the first three terms
  x<-(a[2]/d)
  if(a[3] != 0){
    x <- x +a[3]*((f/d)+l)
  }

  #for the fourth term
  if(t>3){
    x<-x+a[4]
  }
  #initalize some counting variables
  e<-1
  o<-1

  #for all other terms greater than 4
if(t>4){
  for(i in 5:t){

    if(i %% 2 != 0 ){#iff odd
      x<-x+((o+1)*a[i]*f^o)
      o<-o+1

    }
    if(i %% 2 == 0 ){#iff even
      x<-x+a[i]*(((f^(e+1))/d)+(e+1)*(f^e)*l)
      e<-e+1
    }
  }
}
  #some change of variables here
  x<-(x^(-1))

  if(boundedness!='u'){
   M<-pdfQuantileMetalog(a,y,t,bounds=bounds,boundedness='u')
  }
  if(boundedness=='sl'){
    x<-x*exp(-M)
  }
  if(boundedness=='su'){
    x<-x*exp(M)
  }
  if(boundedness=='b'){
    x<-(x*(1+exp(M))^2)/((bounds[2]-bounds[1])*exp(M))
  }

  return(x)
}

#quantile function
quantileMetalog<-function(a,y,t,bounds=c(),boundedness='u'){
  #error check that a is a numeric vector, y is a number between 0,1 and t is greater than a
  #some values for calculation
  f<-(y-0.5)
  l<-log(y/(1-y))

  #for the first three terms
  x<- a[1]+a[2]*l+a[3]*f*l

  #for the fourth term
  if(t>3){
    x<-x+a[4]*f
  }
  #some tarcking variables
  o<-2
  e<-2

  #for all other terms greater than 4
if(t>4){
  for(i in 5:t){
    if(i %% 2 == 0){
        x<-x+a[i]*f^e*l
        e<-e+1
    }
    if(i %% 2 != 0){
       x<-x+a[i]*f^o
       o<-o+1
    }
  }
}
  if(boundedness=='sl'){

    x<-bounds[1]+exp(x)
  }
  if(boundedness=='su'){

    x<-bounds[2]-exp(-x)
  }
  if(boundedness=='b'){

    x<-(bounds[1]+bounds[2]*exp(x))/(1+exp(x))
  }
 return(x)
}

#pdf validation function
#call this feasibility
pdfMetalogValidation <- function(x){
  y<-min(x)
  if(y>=0){
    return('yes')
  }
  if(y<0){
    return('no')
  }

}
