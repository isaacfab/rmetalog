#working on improving the a vector estimation
library(lpSolve)

aVectorsMetalogOLS<-function(x,term_limit,diff){

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

aVectorsMetalogLP<-function(x,term_limit,diff_error=.001,diff_step=0.001){

  A<-data.frame()
  for (i in 1:(term_limit-1)){
    Y<-as.matrix(x[,4:(i+4)])
    Y_neg=-(Y)
    new_Y<-matrix(c(Y[,1],Y_neg[,1]),ncol=2)
    for(c in 2:length(Y[1,])){
      new_Y<-cbind(new_Y,Y[,c])
      new_Y<-cbind(new_Y,Y_neg[,c])
    }
    z<-as.matrix(x$z)
    a<-paste0('a',(i+1))

    #building the constraint matrix
    error_mat<-c()
    for(j in 1:length(Y[,1])){
      front_zeros<-rep(0,(2*(j-1)))
      ones<-c(1,-1)
      trail_zeroes<-rep(0,(2*(length(Y[,1])-j)))
      if(j==1){
        error_vars<-c(ones,trail_zeroes)
      }
      else if(j!=1){
        error_vars<-c(front_zeros,ones,trail_zeroes)
      }
      error_mat<-rbind(error_mat,error_vars)
    }

    new<-cbind(error_mat,new_Y)

    diff_mat<-diffMatMetalog(i+1,diff_step)
    diff_zeros<-c()
    for(t in 1:length(diff_mat[,1])){
      zeros_temp<-rep(0,(2*(length(Y[,1]))))
      diff_zeros<-rbind(zeros_temp,diff_zeros)
    }

    diff_mat<-cbind(diff_zeros,diff_mat)
    #combine the total constraint matrix
    lp_mat<-rbind(new,diff_mat)


    #objective function coeficients
    f.obj<-c(rep(1,(2*length(Y[,1]))),rep(0,(2*(i+1))))
    #constraint matrix
    f.con<-lp_mat
    #symbol vector
    f.dir<-c(rep("==",length(Y[,1])),rep(">=",length(diff_mat[,1])))
    #right hand size for constraints
    f.rhs<-c(z,rep(0,length(diff_mat[,1])))

    #solving the linear program
    lp_sol<-lpSolve::lp("min",f.obj,f.con,f.dir,f.rhs)

    #consolidating solution back into the a vector
    tempLP<-lp_sol$solution[(2*length(Y[,1])+1):(length(lp_sol$solution))]
    temp<-c()
    for(r in 1:(length(tempLP)/2)){
      temp[r]<-tempLP[(r*2-1)]-tempLP[(2*r)]
    }

    #append zeros for term limit
    temp<-c(temp,rep(0,(term_limit-(i+1))))
    if(length(A)==0){
      A<-data.frame(a2=temp)
    }
    if(length(A)!=0){
      A[`a`]<-temp
    }
  }#close the for loop
  return(A)
}#close the function




