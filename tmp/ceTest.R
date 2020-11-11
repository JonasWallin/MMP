##### Test of ce function
## 2020-11-11

library(MMP)
data("sherifdat") # this data is not updated? it is now!

### CEM

# # model on page 536 (table 22.3) in Lang bookchapter (2018)

sherif1 <- ce(y ~ 1 + time, 
                ~ 1 | person, 
                ~ 1 + time | group, 
                emergence = ~ 1 + time, # for CEM we need to include an intercept to get an estimate of sigma
                method = "CEM", 
                data = subset(sherifdat, time >= 0 & time <= 2))

summary.ce(sherif1)

## using old/incorrect data

sherifdat_old=expand.grid(person = 1:3,time = 0:2,
                          group = 1:8)

sherifdat_old$y<-c(4.0,2.7,2.0,3.8,2.6,2.1,2.0,2.4,2.3,1.3,
                   1.5,1.1,2.0,1.6,1.5,1.8,1.5,1.4,2.9,1.9,1.9,2.9,2.1,
                   2.0,1.9,1.7,1.9,3.5,3.4,3.6,4.7,3.8,4.6,4.4,4.4,4.6,
                   4.3,4.3,3.6,4.1,4.0,3.9,4.0,3.9,3.5,4.5,4.0,4.0,1.9,
                   1.8,1.6,1.4,1.6,1.5,2.5,2.5,2.2,2.2,1.6,1.8,
                   
                   1.9,1.8,1.8,3.1,2.6,1.5,5.8,5.2,5.2,4.4,4.4,4.4)

sherifdat_old$condition<-rep(c(1,0),each=36)

sherifdat_old$indv <- 0

indv <- unique(sherifdat_old[,c('person','group')])

for(i in 1:dim(indv)[1]){
  index <- (sherifdat_old$person%in%indv[i,1] )*(sherifdat_old$group%in%indv[i,2] )
  sherifdat_old$indv[index==1] <- i
}

sherif2 <- ce(y ~ 1 + time, 
              ~ 1 | person, 
              ~ 1 + time | group, 
              emergence = ~ 1 + time, # for CEM we need to include an intercept to get an estimate of sigma
              method = "CEM", 
              data = sherifdat_old)

summary.ce(sherif2)



# model with interaction effect on page 536 (table 22.5) in Lang bookchapter (2018) 
sherif3 <- ce(y ~ 1 + time * condition, 
                ~ 1 | person, 
                ~ 1 + time | group, 
                emergence = ~ 1 + time * condition, 
                method = "CEM", 
                data = subset(sherifdat, time >= 0 & time <= 2))

summary.ce(sherif3)


### CEI

sherif4 <- ce(y ~ 1 + time, 
              ~ 1 | person, 
              ~ 1 + time | group, 
              emergence = ~ -1 + time, # for CEI we need to exclude the intercept 
                                       # (why? we still get two indv parameters out?)
                                       # something to do with WLT = F
              method = "CEI", 
              data = sherifdat)

summary.ce(sherif4)
