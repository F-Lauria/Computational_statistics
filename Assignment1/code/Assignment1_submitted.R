#Exercise 1.1.

x = 6
k = 26
m = 1734
n = 203

phyper (x, k, m-k, n, lower.tail = FALSE)

#Exercise 1.2.

incidents<-function(S,n){
  # S is the number of experiments
  # n is the number of the sample size
  i<-0
  out_vec<-rep(i,S)
  lambda<-26/1734
  for(i in 1:S){
    incidents<-rpois(n,lambda)
    number_incidents<-sum(incidents)
    if(number_incidents>=7){
      out_vec[i]=1
    }
  }
  solution<-(sum(out_vec)/S)
  return(solution)
}

set.seed(123)
incidents(10000,203)

#Exercise 1.3.

incidents2<-function(S,n){
  # S is the number of experiments
  # n is the number of the sample size
  i<-0
  out_vec<-rep(i,S)
  lambda<- 26/1734
  for(i in 1:S){
    rateexp=rexp(1,1/lambda)
    incidents<-rpois(n,rateexp)
    number_incidents<-sum(incidents)
    if(number_incidents>=7){
      out_vec[i]=1
    }
  }
  solution<-(sum(out_vec)/S)
  return(solution)
}

set.seed(123)
incidents2(10000,203)



#Exercise 2.1.a

n1<-10
n2<-100
sigma_squared1<-2
sigma_squared2<-1

sp_num<-((n1-1)*sigma_squared1)+((n2-1)*sigma_squared2)
sp_den<-n1+n2-2
SP<-sqrt(sp_num/sp_den)
SE_student<-SP*(sqrt((1/n1)+(1/n2)))

SE_welch<-sqrt((sigma_squared1/n1)+(sigma_squared2/n2))

#Exercise 2.1.c

set.seed(123)
##inputs
sigma_squared1<-2; n1<-c(10,100,200); sigma_squared2<-c(1,2,10);n2<-100;mu1=0;mu2=1;S=10000

#final matrix to put the properties values of all the 9 conditions
all_results<-data.frame(row.names = c("Bias_student","Bias_welch","Variance_student","Variance_welch","MSE_student","MSE_welch","RE_student/welch"))


MonteCarlo=function(S,n1,n2,sigma_squared1,sigma_squared2,mu1,mu2){
  #matrix for calculations
  result<-matrix(nrow=S,ncol=2)
  colnames(result)<-c("SE_student","SE_welch")
  
  for (s in 1:S){
    #computing new variances
    out1<-rnorm(n1,mu1,sqrt(sigma_squared1))
    out2<-rnorm(n2,mu2,sqrt(sigma_squared2))
    out1_var<-var(out1)
    out2_var<-var(out2)
    
    ##SE_student
    sp_num<-((n1-1)*out1_var)+((n2-1)*out2_var)
    sp_den<-n1+n2-2
    SP<-sqrt(sp_num/sp_den)
    SE_student<-SP*(sqrt((1/n1)+(1/n2)))
    SE_student<-SP*(sqrt((1/n1)+(1/n2)))
    ##SE_Welch
    SE_welch<-sqrt((out1_var/n1)+(out2_var/n2))
    
    ##Storing results
    result[s,]<-c(SE_student,SE_welch)
  }
  
  return(result)
}


Evaluation=function(result,n1,n2,sigma_squared1,sigma_squared2){
  #True value parameter
  teta<-sqrt((sigma_squared1/n1)+(sigma_squared2/n2))
  
  #Bias
  MC_mean<-apply(result,2,mean)
  Bias_student<-teta-MC_mean[1]
  Bias_welch<-teta-MC_mean[2]
  
  #Variance
  MC_var<-apply(result,2,var) 
  Variance_Student<-MC_var[1]
  Variance_Welch<-MC_var[2]
  
  #MSE
  MC_MSE_student<-(Bias_student^2)+MC_var[1]
  MC_MSE_welch<-(Bias_welch^2)+MC_var[2]
  
  #Relative efficiency
  MC_RE<-MC_MSE_student/MC_MSE_welch
  
  #Storing Values
  resultdf<-data.frame(Bias_student,Bias_welch,Variance_Student,Variance_Welch,MC_MSE_student,MC_MSE_welch,MC_RE,row.names = "case 0" )
  return(resultdf)
}



#Create the final table thanks to the functions above 
for (i in 1:3){
  for (j in 1:3){
    MC_result<-MonteCarlo(S,n1[i],n2,sigma_squared1,sigma_squared2[j],mu1,mu2)
    appendThat<-Evaluation(MC_result,n1[i],n2,sigma_squared1,sigma_squared2[j])
    all_results<-rbind.data.frame(all_results,appendThat)
    
  }
}  
    

