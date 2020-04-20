#Load data
VetSuicides_2005 <- read.csv("data/VetSuicides_2005.csv")
VetSuicides_2011 <- read.csv("data/VetSuicides_2011.csv")
## store the interested variables
vet_rate_2005<-VetSuicides_2005$vet_rate 
vet_rate_2011<-VetSuicides_2011$vet_rate
#--------------------------------- Q1.a ---------------------------------------- #Checking distribution
{qqnorm(vet_rate_2005,main = "Q-Q plot 2005 Vet_rate")
  qqline(vet_rate_2005)}
{qqnorm(vet_rate_2011,main = "Q-Q plot 2011 Vet_rate") 
  qqline(vet_rate_2011)}
#Computing means
sample_mean_2011 <- mean(vet_rate_2011) 
sample_mean_2005 <- mean(vet_rate_2005)
#Dataset to compute
Diff <- vet_rate_2011 - vet_rate_2005
#Get a histogram of the difference in means.
hist(Diff, xlim = c(-30,30)) 
abline(v = mean(Diff), col=2, lwd=2)
#--------------------------------- Q1.b ---------------------------------------- #Performing t-test
ttest <- t.test(vet_rate_2011,vet_rate_2005, var.equal = TRUE, alternative="greater",
                paired=TRUE)
ttest
#--------------------------------- Q1.c ----------------------------------------
## NON PRAMETRIC BOOTSTRAP:
non_parametric_simulation<-function(B,vet_rate_2011,vet_rate_2005){
  result<-matrix(nrow=B,ncol=2) 
  colnames(result)<-c("Mean_2011","Mean_2005")
for (i in 1:B){
  indices<-sample(c(1:50),replace=TRUE)
  bootsample_2011<-vet_rate_2011[indices]
  mean_2011 <- mean(bootsample_2011)
  bootsample_2005<-vet_rate_2005[indices]
  mean_2005 <- mean(bootsample_2005)
  result[i,]<-c(mean_2011,mean_2005)
}
  return(result) 
}

set.seed(12345) 
result<-non_parametric_simulation(15000,vet_rate_2011,vet_rate_2005)

# sampling distribution of the paired different in means
bootstrap_distribution<-(result[,1]-result[,2])

#Plotting bootstrap distribution
hist(bootstrap_distribution,main = "Non-parametric bootstrap distribution", xlab = "Mean difference")

#Computing Bias
obs_mean<-mean(vet_rate_2011)-mean(vet_rate_2005) 
boot_mean<-mean(bootstrap_distribution) 
bias<-boot_mean-obs_mean # NO BIAS

#Confidence intervals
#Get the lenght
N<-length(vet_rate_2005) 
M<-length(vet_rate_2011)

SE_bootstrap<-sd(bootstrap_distribution) 
alfa <- 0.05
tcv <- qt(1-alfa,N+M-2)

upper <- boot_mean+(tcv*SE_bootstrap)
upper_quant <- quantile(bootstrap_distribution,0.95)

#--------------------------------- Q1.d ----------------------------------------
#Permutation
all_data<-matrix(nrow = 50, ncol = 2)
all_data[1:50,1]<-vet_rate_2011
all_data[1:50,2]<-vet_rate_2005
colnames(all_data)<-c("vet_rate_2011","vet_rate_2005")

set.seed(12345)

# Permutation procedure:
B <- 9999 #nr of permutation samples

out <- matrix(nrow = B, ncol = 2) 
groups<-matrix(nrow = 50, ncol = 2)


colnames(groups)<-c("group_1","group_2")

for (i in 1:B){
  for(row in 1:50){ 
    groups[row,]<-sample(all_data[row,],2)
  }

  mean_group1<-mean(groups[,1])
  mean_group2<-mean(groups[,2])
  out[i,] <- c(mean_group1,mean_group2)
}

#Calculate p-value
diff<-out[,1]-out[,2]
obsdiff<-mean(vet_rate_2011)-mean(vet_rate_2005)
pval<-(1+sum(diff>(obsdiff)))/(B+1)

#Creating the plot.
a <- hist(diff,50)


#Permutation distribution
barcols = a$breaks
barcols[a$breaks<round(obsdiff,0)] = 0
barcols[a$breaks>=round(obsdiff,0)] = 2
hist(diff,50, xlab = "",
    main = "Permutation distribution of the mean difference (2011-2005)",
    col=barcols,xlim = c(-6,6))
abline(v=obsdiff, col =2, lwd = 2)

#--------------------------------- Q2.2 ---------------------------------------- 
#Monte carlo simulation

#inputs
sigma_squared2 = c(1, 2, 10) ; n1 = c(10, 100, 200) 
n2<- 100 ; sigma_squared1 = 2
S = 1000; B = 99
set.seed(3)


MonteCarlo=function(S,B,n1,n2,sigma_squared1,sigma_squared2,mu1,mu2){

  result <- matrix(nrow = S ,ncol = 3 )
  colnames(result)<-c("student p", "welch p", "perm")
  for (s in 1:S){
    #Getting 2 sample data.
    out1<-matrix(rnorm(n1,mu1,sqrt(sigma_squared1)))
    out2<-matrix(rnorm(n2,mu2,sqrt(sigma_squared2)))
    
    #Getting the Welch and Student
    Welcht <- t.test(out1, out2, alternative = "two.sided")
    Studentt <- t.test(out1, out2, var.equal = T, alternative = "two.sided")
  
    ## Permutation
  
    #Getting a combined data set.
    all_out <- rbind(out1,out2)
    #Getting the difference in means.
    obs_diff <- mean(out1) - mean(out2)
  
    # Permutation procedure:
    out <- matrix(nrow = B, ncol = 2)
  
    for (i in 1:B){
      permnrs <- sample(n1+n2) #this is the re-sampling without replacement
      xbar_group1 <- mean(all_out[permnrs[1:n1],1]) # in this way we label the sample
      xbar_group2 <- mean(all_out[permnrs[(n1+1):(n1+n2)],1])
      out[i,] <- c(xbar_group1,xbar_group2)
    }
  #calculate p-value
    diff<-out[,1]-out[,2]
    pval<-(1 + sum(abs(diff) > abs(obs_diff)))/(B+1) # larger or different
  
  ##Storing results
  
  result[s,]<-c(Studentt$p.value, Welcht$p.value, pval)
  }
return(result) 
}


#Null hypothesis, Computing Levels, 9 cases.

mu1=0;mu2=0

Null_Hypothesis_level<- data.frame(row.names = c("n1","sigma_sq_2","Student lvl",
                                                 "Welch LvL", "Permutation LvL"))

for (i in 1:3){
  for (j in 1:3){
    result <- MonteCarlo(S,B,n1[i],n2,sigma_squared1,sigma_squared2[j],mu1,mu2)
    
    #Getting levels
    levelS <- sum(result[,1]<0.05)/S
    levelW <- sum(result[,2]<0.05)/S
    levelP <- sum(result[,3]<0.05)/S
    
    Result_level <- data.frame(n1[i],sigma_squared2[j],levelS,levelW,levelP)
    
    Null_Hypothesis_level <- rbind.data.frame(Null_Hypothesis_level,
                                              Result_level)
  }
}

#Alternative hypothesis, computing Power, 9 cases.
mu1=0;mu2=1

Alt_Hypothesis_power<- data.frame(row.names = c("n1","sigma_sq_2","Student Power",
                                                "Welch Power", "Permutation Power"))

for (i in 1:3){
  for (j in 1:3){
    result <- MonteCarlo(S,B,n1[i],n2,sigma_squared1,sigma_squared2[j],mu1,mu2)
    #Getting power
    PowerS <- sum(result[,1]<0.05)/S
    PowerW <- sum(result[,2]<0.05)/S
    PowerP <- sum(result[,3]<0.05)/S
    
    Result_level <- data.frame(n1[i],sigma_squared2[j],PowerS, PowerW, PowerP)
    Alt_Hypothesis_power <- rbind.data.frame(Alt_Hypothesis_power,
                                             Result_level)
  }
}
  
  
  
#Matrix demonstration Level
Null_Hypothesis_level

#Matrix demonstration Power
Alt_Hypothesis_power
  

