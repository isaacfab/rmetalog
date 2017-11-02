#will use this file for inital metalog functions
#example data for now

#error checking here, valid values for the inputs
#if prob data is input need to check it for validity
#boundedness
#must be n,sl,su or b
#must have a numeric upper bound
#must have a numeric lower bound

#bounds
#lower and upper bounds must be less/greater than all x values

#prob
#each entry must have a quantile
#all entries need to be numeric between 0 and 1
#data must be contiguous

#x
#each entry needs to be numeric
#data must be contiguious
#length of data must be at least 2

#term_limit
#must be an integer in the range of 5 to n
#must be less than the length of x
rMetalog <- function(x,step_len=.01,probs=0,term_limit=16,bounds=c(),boundedness='u') {

#create a list to hold all the objects
myList<-list()
#################inital error checking################
if(class(x)!='numeric'){
  return(print('Error: input x must be a numeric vector!'))
}

###############handle the probabilites###############
#this also converts x as a data frame
   if(probs == 0){
     x<-MLprobs(x)
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
  x$y1<-1
  x$y2<-(log(x$probs/(1-x$probs)))
  x$y3<-(x$probs-0.5)*x$y2
  x$y4<-x$probs-0.5

#####complete the values through the term limit#####
for (i in 5:(term_limit)){

    y<-paste0('y',i)
    if(i %% 2 != 0){
     x[`y`]<-x$y4^(i%/%2)
    }
    if(i %% 2 == 0){
     z<-paste0('y',(i-1))
     x[`y`]<-x$y2*x[`z`]
    }
}
  myList$Y<-x

###########build a vectors for each term###########
  A<-data.frame()
  for (i in 1:(term_limit-1)){
    Y<-as.matrix(x[,4:(i+4)])
    z<-as.matrix(x$z)
    a<-paste0('a',(i+1))
    #add error catching here for non invertable
    temp<-((solve(t(Y)%*%Y) %*% t(Y)) %*% z)
    temp<-c(temp,rep(0,(term_limit-(i+1))))
    if(length(A)==0){
      A<-data.frame(a2=temp)
    }
    if(length(A)!=0){
      A[`a`]<-temp
    }
  }
  myList$A<-A

##############build the metalog m and M dataframes###############
y<-seq(step_len,(1-step_len),step_len)

Mh<-data.frame()
for(i in 2:term_limit){
  a_name<-paste0('a',i)
  m_name<-paste0('m',i)
  M_name<-paste0('M',i)

#build pdf
  m<-pdfMetalog(myList$A[`a_name`][,1],y[1],term_limit,bounds=bounds,boundedness=boundedness)

  for(j in 2:length(y)){
    temp<-pdfMetalog(myList$A[`a_name`][,1],y[j],term_limit,bounds=bounds,boundedness=boundedness)
    m<-c(m,temp)
  }

#build inverse
  M<-pdfQuantileMetalog(myList$A[`a_name`][,1],y[1],term_limit,bounds=bounds,boundedness=boundedness)

  for(j in 2:length(y)){
    temp<-pdfQuantileMetalog(myList$A[`a_name`][,1],y[j],term_limit,bounds=bounds,boundedness=boundedness)
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
}


myList$M<-Mh

#########pdf validation################
y<-c()
for(i in 2:term_limit){
  m<-paste0('m',i)
  x<-pdfMetalogValidation(myList$M[`m`])
  y<-c(y,x)

}

myList$Validation<-y


return(myList)
}

myMetalog <- rMetalog(x$FishSize,step_len = .001,bounds=c(0,60),boundedness = 'u',term_limit = 16)
#myMetalog <- rMetalog(data.frame(sort(myInterestingData$time_diff_num)),step_len = .001,boundedness = 'sl')


#need a fuction that retuns a cdf function from an emperical input


#will us this to build the base plots
#myfun1 <- function(x){
#  x + 2.5
#}

#myfun2 <- function(x){
#  -0.5*x + 4.5
#}

#x<-c(0,0,0.5,3,4,4)
#y<-c(0,2.5,3,3,2.5,0)

#curve(myfun1(x),from=-2,to=5,col="blue",xlab="x1",ylab="x2")
#abline(h=0, v=0, col = "gray60",lwd=5)
#curve(myfun2(x),add=TRUE,from=-1,to=10,col="red")
#abline(v=4, col = "red")
#abline(h=3, col = "red")
#polygon(x,y,col="blue")
