#working on improving the a vector estimation
library(lpSolve)

#this function will coduct a two stage metalog optimization. It will attempt an OLS original method and if it fails or is not valid will run the LP method
a_vector_OLS_and_LP<-function(myList,term_limit,term_lower_bound,bounds,boundedness,diff_error=.001,diff_step=0.001){
  #some place holder values
  A<-data.frame()
  c_a_names<-c()
  c_m_names<-c()
  Mh<-data.frame()
  Validation<-data.frame()
  #cat('Building the metalog distributions now')
  #pb<-progress::progress_bar$new(total=(term_limit-(term_lower_bound-1)))
  for (i in term_lower_bound:term_limit){
    #pb$tick()
    #fix this matrix with a Y object
    Y<-as.matrix(myList$Y[,1:(i)])
    z<-as.matrix(myList$dataValues$z)
    y<-myList$dataValues$probs

    a<-paste0('a',(i))
    m_name<-paste0('m',i)
    M_name<-paste0('M',i)
    c_m_names<-c(c_m_names,m_name,M_name)
    c_a_names<-c(c_a_names,a)
    #try to use the OLS approach
    temp<-try(((solve(t(Y)%*%Y) %*% t(Y)) %*% z),silent = TRUE)
    #if temp is not a valid numeric vector then move to LP
    if(class(temp)!='matrix'){
        temp<-a_vector_LP(myList,
                          term_limit=i,
                          term_lower_bound=i,
                          diff_error=diff_error,
                          diff_step=diff_step)
    }

    temp<-c(temp,rep(0,(term_limit-(i))))

    #get the list and quantile values back for validation
    tempList<-pdf_quantile_builder(temp,
                                   y,
                                   term_limit=i,
                                   bounds=bounds,
                                   boundedness=boundedness)

    #if it not a valid pdf run and the OLS version was used the LP version
    if(tempList$valid=='no' & class(temp)=='numeric'){
          temp<-a_vector_LP(myList,
                            term_limit=i,
                            term_lower_bound=i,
                            diff_error=diff_error,
                            diff_step=diff_step)
          temp<-c(temp,rep(0,(term_limit-(i))))
          #get the list and quantile values back for validation
          tempList<-pdf_quantile_builder(temp,
                                         y,
                                         term_limit=i,
                                         bounds=bounds,
                                         boundedness=boundedness)
    }

    if(length(Mh)!=0){
      Mh<-cbind(Mh,tempList$m)
      Mh<-cbind(Mh,tempList$M)
    }
    if(length(Mh)==0){
      Mh<-as.data.frame(tempList$m)
      Mh<-cbind(Mh,tempList$M)
    }
    if(length(A)!=0){
      A<-cbind(A,temp)
    }
    if(length(A)==0){
      A<-as.data.frame(temp)
    }
    tempValidation<-data.frame(term=i,valid=tempList$valid)
    Validation<-rbind(Validation,tempValidation)
  }#close the for loop
  colnames(A)<-c_a_names
  colnames(Mh)<-c_m_names
  myList$A<-A
  myList$M<-Mh
  myList$M$y<-tempList$y
  myList$Validation<-Validation
  return(myList)
}

a_vector_LP<-function(myList,term_limit,term_lower_bound,diff_error=.001,diff_step=0.001){

  A<-data.frame()
  cnames<-c()

  for (i in term_lower_bound:term_limit){


    Y<-as.matrix(myList$Y[,1:(i)])
    z<-as.matrix(myList$dataValues$z)

    #bulding the objective function using abs value LP formulation
    Y_neg=-(Y)
    new_Y<-matrix(c(Y[,1],Y_neg[,1]),ncol=2)
    for(c in 2:length(Y[1,])){
      new_Y<-cbind(new_Y,Y[,c])
      new_Y<-cbind(new_Y,Y_neg[,c])
    }
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
    #temp<-c(temp,rep(0,(term_limit-(i))))

    if(length(A)!=0){
      A<-cbind(A,temp)
    }
    if(length(A)==0){
      A<-as.data.frame(temp)
    }

  }#close the for loop
  colnames(A)<-cnames
  return(temp)
}#close the function




