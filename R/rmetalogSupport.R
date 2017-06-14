#Supporting functions called inside the metalog function call

#build the quantiles through a base function
MLprobs <- function(x) {
  l <- length(x[,1])
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

pdfMetalog<-function(a,y,t){
  #error check that a is a numeric vector, y is a number between 0,1 and t is greater than a
  #some values for calculation
  d<-y*(1-y)
  f<-(y-0.5)
  l<-log(y/(1-y))

  #initiate pdf

  x<-(a[2]/d)+a[4]
  if(a[3] != 0){
    x <- x +a[3]*((f/d)+l)
  }

  #initalize some counting variables
  e<-1
  o<-1
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
  return((x^(-1)))
}

#Inverse lookup
#cdf
pdfInvMetalog<-function(a,y,t){
  #error check that a is a numeric vector, y is a number between 0,1 and t is greater than a
  #some values for calculation
  f<-(y-0.5)
  l<-log(y/(1-y))

  x<- a[1]+a[2]*l+a[3]*f*l+a[4]*f

  #some tarcking variables
  o<-2
  e<-2

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
 return(x)
}

#pdf validation function
#call this feasibility
pdfMetalogValidation <- function(x){
  y<-min(x)
  if(y>0){
    return('yes')
  }
  if(y<0){
    return('no')
  }

}
