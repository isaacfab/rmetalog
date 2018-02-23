#working on improving the a vector estimation
library(lpSolve)

aVectorMetalogOLS<-function(x,term_limit){

  A<-data.frame()
  for (i in 1:(term_limit-1)){
    Y<-as.matrix(x[,4:(i+4)])
    z<-as.matrix(x$z)
    a<-paste0('a',(i+1))
    #add error catching here for non invertable
    tempOG<-((solve(t(Y)%*%Y) %*% t(Y)) %*% z)
    temp<-c(temp,rep(0,(term_limit-(i+1))))
    if(length(A)==0){
      A<-data.frame(a2=temp)
    }
    if(length(A)!=0){
      A[`a`]<-temp
    }
  }
return(A)
}

aVectorMetalogLP<-function(x,term_limit){

  A<-data.frame()
  for (i in 1:(term_limit-1)){
    Y<-as.matrix(x[,4:(i+4)])
    z<-as.matrix(x$z)
    #objective function coeficients
    f.obj<-0
    #constraint matrix
    f.con<-0
    #symbol vector
    f.dir<-c()
    #right hand size for constraints
    f.rhs<-0

    #lp_sol<-lpSolve::lp("min",f.obj,f.con,f.dir,f.rhs,all.bin=TRUE)
  }#close the for loop
  return(A)
}#close the function



#temp<-aVectorMetalog(x,20)
