##### R libraries required
library("truncnorm")
# library("doParallel") ####if you can run in parallel (i.e. your machine has enough cores)
# library("erer") ### for writing output at the end

#### First read in empirical dataset which has parameter values
paramdis = read.csv(choose.files(),stringsAsFactors = FALSE)

head(paramdis)
str(paramdis)

#cht = proportional change in impact
#chc = proportional change in control
#tcd = proportional difference between impact and control (impact-control) in before period

paramdis = data.frame(paramdis)
paramdis$Chc = as.numeric(paramdis$Chc)
paramdis$Cht = as.numeric(paramdis$Cht)
paramdis$tcd= as.numeric(paramdis$tcd )

runs=1000 ### number of unique runs in simulation
multb1 = 1000 #### number of repetitions of each run

### find 25% and 75% quantiles for each parameter
quantschc = quantile(paramdis$Chc,na.rm=TRUE)[c(2,4)]
quantscht =quantile(paramdis$Cht,na.rm=TRUE)[c(2,4)]
quantstcd =quantile(paramdis$tcd,na.rm=TRUE)[c(2,4)]

### restrict parameter values to InterQuartile Range
Chc2 = paramdis$Chc[which(paramdis$Chc>=quantschc[1] & paramdis$Chc<=quantschc[2])] 
Cht2 = paramdis$Cht[which(paramdis$Cht>=quantscht[1] & paramdis$Cht<=quantscht[2])]
tcd2 = paramdis$tcd[which(paramdis$tcd>=quantstcd[1] & paramdis$tcd<=quantstcd[2])] 

#### sample from parameter values 1000 times with replacement to generate 1000 unique parameter runs
#### then repeat these 1000 times (1000 repetitions)
cht = rep(round(sample(Cht2[which(is.na(Cht2)==FALSE)],runs,replace=TRUE),digits=3),multb1)
chc = rep(round(sample(Chc2[which(is.na(Chc2)==FALSE)],runs,replace=TRUE),digits=3),multb1)
tcd = rep(round(sample(tcd2[which(is.na(tcd2)==FALSE)],runs,replace=TRUE),digits=3),multb1)

### sample from truncated normal distribution to get values of between site coefficient of variation (CV) for each of the 1000 runs
svar1= rep(rtruncnorm(runs,0,0.2,0.1,0.05),multb1)

### put parameters together in dataframe ready for simulation
params = data.frame(cbind(chc,cht,tcd,svar1))

write.csv(params,"params.csv",row.names=FALSE) ### for reference

#### function for simulation
###zeta = repetition number

supersimsc <- function(zeta){
  
  #### input parameters from dataframe we created
  svar=params$svar1[zeta]
  
  ### set lambda to 50 for population abundance
  startvalt = 50
  
  ### set values for change in impact, control and initial difference between impact and control
  cht1 = (params$cht[zeta]-1)*startvalt
  chc1 = (params$chc[zeta]-1)*startvalt
  tcd1 = (params$tcd[zeta]-1)*startvalt
  
  
  ### Possible combinations of sites and time steps
  sitesvals1t = c(2,5,10,25,50) ###possble numbers of impact sites sampled
  sitesvals1c = c(2,5,10,25,50) ###possble numbers of control sites sampled
  yrst = c(2,4,6,8,10) ###possble numbers of time steps sampled and simulated
  ##### beware amount of memory this will require, data produced = ~9GB!
  
  ### create a matrix to capture effect sizes for each design and all possible combinations of spatial replication (5x5) and time steps sampled (x4)
  results = matrix(ncol=14,nrow=length(sitesvals1t)*length(sitesvals1c)*length(c(2,4,6,8,10)))
  z=1
  
  for (tt1 in 1:length(yrst)){ ###loop through all possible numbers of time steps simulated
    
    ###set number of time steps to simulate in each period 
    ###(i.e. if maxyears = 10, 10 before + 10 after = 20 time steps total)
    maxyears = yrst[tt1]
    
    ### first generate time steps for impact and control group in before period
    ### then generate same number of time steps, but with a lambda adjusted by cht1 and chc1 (change in impact and control) 
    year1tat = c(rpois(maxyears,startvalt),rpois(maxyears,startvalt+cht1))
    year1tac = c(rpois(maxyears,startvalt),rpois(maxyears,startvalt+chc1))
    
    ### loop through all these combinations of sites and time steps sampled and generate effect sizes
    
    for(st in 1:length(sitesvals1t)){
      for(sc in 1:length(sitesvals1c)){
        
        s1=sitesvals1t[st] # no. impact sites sampled
        s2=sitesvals1c[sc] # no. control sites sampled
        tt=yrst[tt1] # no. time steps sampled
        
        #### create matrices to capture site data for each time step and group
        
        ### randomised controlled trial design 
        temptranci = matrix(ncol=tt,nrow=s1)### After Impact
        tempcranci = matrix(ncol=tt,nrow=s2) ### After Control
        
        ### non-randomised designs: baci (using all 4 matrices) or BA (using only first 2 matrices)
        ### or CI (using 2nd and 4th matrices) or After (using only 2nd matrix)
        temptb = matrix(ncol=tt,nrow=s1) ### Before Impact
        tempta = matrix(ncol=tt,nrow=s1) ### After Impact
        tempcb = matrix(ncol=tt,nrow=s2) ### Before Control
        tempca = matrix(ncol=tt,nrow=s2) ### After Control
        
        for(t1 in 1:tt){ ### for every time step, sample sites
          
          #### use pmax to prevent negative values
          #### generate sites using normal distribution with a mean taken from the time step value
          #### svar is the coefficient of variation (CV) we parameterised earlier and is multiplied by the mean value to give the sd of the normal distribution
          
          ### randomised controlled trial design 
          temptranci[,t1] = pmax((rnorm(s1, year1tat[maxyears+t1],year1tat[maxyears+t1]*svar)),0) #After impact
          tempcranci[,t1] = pmax((rnorm(s2, year1tac[maxyears+t1],year1tac[maxyears+t1]*svar)),0) #After control
          
          ### After or CI or Before-After or BACI design data
          ###note addition of tcd1 for control data to reflect non-random site selection so 
          ###impact and control groups differ on average by tcd1
          
          tempta[,t1] = pmax((rnorm(s1, year1tat[maxyears+t1],year1tat[maxyears+t1]*svar)),0) #After impact
          temptb[,t1] = pmax((rnorm(s1, year1tat[maxyears-t1+1],year1tat[maxyears-t1+1]*svar)),0) #Before impact
          
          tempca[,t1] = pmax((rnorm(s2, pmax(year1tac[maxyears+t1]+tcd1,0),pmax(year1tac[maxyears+t1]+tcd1,0)*svar)),0) #After control
          tempcb[,t1] = pmax((rnorm(s2, pmax(year1tac[maxyears-t1+1]+tcd1,0),pmax(year1tac[maxyears-t1+1]+tcd1,0)*svar)),0) #Before control
          
        }
        
        ### find means for impact and control in each period and use to calculate effect size for each design
        ### store them along with how many time steps and sites were sampled and repeat
        
        results[z,1] <- mean(temptranci) - mean(tempcranci) ###RCT
        results[z,2] <- mean(tempta) - mean(temptb)  ### BA
        results[z,3] <- mean(tempta) - mean(tempca) ### CI
        results[z,4] <- (mean(tempta) -  mean(tempca)) - (mean(temptb)-mean(tempcb))  ###BACI
        results[z,5] <- (mean(tempta[,tt]) -  mean(tempta[,1])) ### After using final time step sampled - first time step sampled
        results[z,6] <- tt ###number of time steps sampled
        results[z,7] <- s1 ###number of impact sites sampled
        results[z,8] <- s2 ###number of control sites sampled
        
        ###true change over simulated number of time steps
        results[z,9] <- (mean(year1tat[(maxyears+1):(maxyears+tt)]) -  mean(year1tat[(maxyears-tt+1):maxyears])) -
          (mean(year1tac[(maxyears+1):(maxyears+tt)])-mean(year1tac[(maxyears-tt+1):maxyears]))
        
        
        results[z,10] <- qnorm(0.975)*sqrt(((((length(c(temptranci))-1)*var(c(temptranci)) + (length(c(tempcranci))-1)*var(c(tempcranci)))/(length(c(tempcranci))+length(c(temptranci))-2))/length(c(tempcranci)))
                                           + ((((length(c(temptranci))-1)*var(c(temptranci)) + (length(c(tempcranci))-1)*var(c(tempcranci)))/(length(c(tempcranci))+length(c(temptranci))-2))/length(c(temptranci))))###RCT Error
        results[z,11] <- qnorm(0.975)*sqrt(((((length(c(tempta))-1)*var(c(tempta)) + (length(c(temptb))-1)*var(c(temptb)))/(length(c(temptb))+length(c(tempta))-2))/length(c(temptb)))
                                           + ((((length(c(tempta))-1)*var(c(tempta)) + (length(c(temptb))-1)*var(c(temptb)))/(length(c(temptb))+length(c(tempta))-2))/length(c(tempta))))### BA Error
        results[z,12] <- qnorm(0.975)*sqrt(((((length(c(tempta))-1)*var(c(tempta)) + (length(c(tempca))-1)*var(c(tempca)))/(length(c(tempca))+length(c(tempta))-2))/length(c(tempca)))
                                           + ((((length(c(tempta))-1)*var(c(tempta)) + (length(c(tempca))-1)*var(c(tempca)))/(length(c(tempca))+length(c(tempta))-2))/length(c(tempta))))### CI error
        results[z,13] <- qnorm(0.975)*sqrt(((((length(c(tempta))-1)*var(c(tempta)) + (length(c(tempca))-1)*var(c(tempca)) + (length(c(temptb))-1)*var(c(temptb))+ (length(c(tempcb))-1)*var(c(tempcb)))/(length(c(tempca))+length(c(tempta))+length(c(tempcb))+length(c(temptb))-4))/length(c(tempca)))
                                           + ((((length(c(tempta))-1)*var(c(tempta)) + (length(c(tempca))-1)*var(c(tempca)) + (length(c(temptb))-1)*var(c(temptb))+ (length(c(tempcb))-1)*var(c(tempcb)))/(length(c(tempca))+length(c(tempta))+length(c(tempcb))+length(c(temptb))-4))/length(c(tempta)))
                                           + ((((length(c(tempta))-1)*var(c(tempta)) + (length(c(tempca))-1)*var(c(tempca)) + (length(c(temptb))-1)*var(c(temptb))+ (length(c(tempcb))-1)*var(c(tempcb)))/(length(c(tempca))+length(c(tempta))+length(c(tempcb))+length(c(temptb))-4))/length(c(tempcb)))
                                           +((((length(c(tempta))-1)*var(c(tempta)) + (length(c(tempca))-1)*var(c(tempca)) + (length(c(temptb))-1)*var(c(temptb))+ (length(c(tempcb))-1)*var(c(tempcb)))/(length(c(tempca))+length(c(tempta))+length(c(tempcb))+length(c(temptb))-4))/length(c(temptb)))) ##BACI error
        results[z,14] <- qnorm(0.975)*sqrt(((((length(c(tempta[,tt]))-1)*var(c(tempta[,tt])) + (length(c(tempta[,1]))-1)*var(c(tempta[,1])))/(length(c(tempta[,1]))+length(c(tempta[,tt]))-2))/length(c(tempta[,1])))
                                           + ((((length(c(tempta[,tt]))-1)*var(c(tempta[,tt])) + (length(c(tempta[,1]))-1)*var(c(tempta[,1])))/(length(c(tempta[,1]))+length(c(tempta[,tt]))-2))/length(c(tempta[,tt]))))### After design error
        
        
        
        
        
        z=z+1
      }
    }
  }
  
  ###package up results
  return(list(cbind(results)))
  
}


### test how long it will take serially (without parallelisation) - test with 100 repetitions and multiply up to get a rough estimate of total time (multiply time by 10000)
samp = sample(1:(nrow(params)),100,replace=FALSE)

system.time({
  sapply(samp,supersimsc)
})


### we ran this in parallel on a windows machine with the following code
ncores = 8 ### put in the number of cores you want to use - on our Hyper-threaded Core i9 (20 core) PC 
###we found the best performance was achieved at 8 cores (probably as true number of cores is 10)

c1 <- makeCluster(ncores) ###create a cluster of cores
clusterExport(c1,c('params')) ### import the empirical dataset that it will need to run simulation
clusterEvalQ(c1,library('truncnorm')) ### import any libraries it may need

y <- (1:nrow(params)) ###find number of repetitions (1x10^6)

system.time(
  devall <- parSapplyLB(c1,y,supersimsc) #### time it with system.time if desired and run in parallel with load balancing (LB)
)

stopCluster(c1) ### make sure to stop cluster at the end
gc()

### parallel processing generates list output, advise writing out in chunks as we generate a lot of data!
nchunks=4
chunks = split(1:nrow(params), 1:nchunks)
str(chunks)

for(i in 1:nchunks){
write.list(devall[chunks[[i]]],file=paste("simresults",i,".csv",sep=""),row.names=FALSE)
}

