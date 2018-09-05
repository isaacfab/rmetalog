#working on improving the a vector estimation
library(lpSolve)

aVectorsMetalogOLS<-function(Ymat,x,term_limit,term_lower_bound,diff){

  A<-data.frame()
  for (i in 1:(term_limit-1)){
    #fix this matrix with a Y object
    Y<-as.matrix(Ymat[,1:(i)])
    z<-as.matrix(x$z)
    a<-paste0('a',(i))
    #add error catching here for non invertable
    tempOG<-((solve(t(Y)%*%Y) %*% t(Y)) %*% z)
    temp<-c(temp,rep(0,(term_limit-(i))))
    ############################
    #this needs to be updated
    ############################
    if(length(A)==0){
      A<-data.frame(a2=temp)
    }
    if(length(A)!=0){
      A[`a`]<-temp
    }
  }
return(A)
}

aVectorsMetalogLP<-function(Ymat,x,term_limit,term_lower_bound,diff_error=.001,diff_step=0.001){

  A<-data.frame()
  cnames<-c()
  print('Building the metalog distributions now', row.names=FALSE)
  pb<-progress::progress_bar$new(total=(term_limit-(term_lower_bound-1)))
  for (i in term_lower_bound:term_limit){
    pb$tick()

    Y<-as.matrix(Ymat[,1:(i)])

    #bulding the objective function using abs value LP formulation
    Y_neg=-(Y)
    new_Y<-matrix(c(Y[,1],Y_neg[,1]),ncol=2)
    for(c in 2:length(Y[1,])){
      new_Y<-cbind(new_Y,Y[,c])
      new_Y<-cbind(new_Y,Y_neg[,c])
    }
    z<-as.matrix(x$z)
    a<-paste0('a',(i))
    cnames<-c(cnames,a)
    #building the constraint matrix
    error_mat<-c()
    for(j in 1:nrow(Y)){
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

    diff_mat<-diffMatMetalog(i,diff_step)
    diff_zeros<-c()
    for(t in 1:length(diff_mat[,1])){
      zeros_temp<-rep(0,(2*(length(Y[,1]))))
      diff_zeros<-rbind(zeros_temp,diff_zeros)
    }

    diff_mat<-cbind(diff_zeros,diff_mat)
    #combine the total constraint matrix
    lp_mat<-rbind(new,diff_mat)


    #objective function coeficients
    f.obj<-c(rep(1,(2*length(Y[,1]))),rep(0,(2*(i))))
    #constraint matrix
    f.con<-lp_mat
    #symbol vector
    f.dir<-c(rep("==",length(Y[,1])),rep(">=",length(diff_mat[,1])))
    #right hand size for constraints
    f.rhs<-c(z,rep(diff_error,length(diff_mat[,1])))

    #solving the linear program
    lp_sol<-lpSolve::lp("min",f.obj,f.con,f.dir,f.rhs)

    #consolidating solution back into the a vector
    tempLP<-lp_sol$solution[(2*length(Y[,1])+1):(length(lp_sol$solution))]
    temp<-c()
    for(r in 1:(length(tempLP)/2)){
      temp[r]<-tempLP[(r*2-1)]-tempLP[(2*r)]
    }

    #append zeros for term limit
    temp<-c(temp,rep(0,(term_limit-(i))))

    if(length(A)!=0){
      A<-cbind(A,temp)
    }
    if(length(A)==0){
      A<-as.data.frame(temp)
    }

  }#close the for loop
  colnames(A)<-cnames
  return(A)
}#close the function




