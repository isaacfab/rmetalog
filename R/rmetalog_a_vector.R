#working on improving the a vecotr estimation
library(lpSolve)

aVectorMetalogOG<-function(x,term_limit){

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

aVectorMetalogOG<-function(x,term_limit){

  A<-data.frame()
  for (i in 1:(term_limit-1)){
    Y<-as.matrix(x[,4:(i+4)])
    z<-as.matrix(x$z)
    #column sums
    y_cs<-colSums(Y)
    z_cs<-sum(z)
    #objective function coeficients
    f.obj<-c(y_cs,-(z_cs))

  }#close the for loop
  return(A)
}#close the function

  lp_sol<-lpSolve::lp("max",f.obj,f.con,f.dir,f.rhs,all.bin=TRUE)

temp<-aVectorMetalog(x,20)
