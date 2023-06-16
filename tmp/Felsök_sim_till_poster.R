library(MMP)

data <- readRDS("tmp/Postervar_delta_6.Rdata")

data6 <- data[[2]]
data6_miss <- data6[-data$indices[[3]],]

m1 <- ce(y ~ 1+time, 
         ~ 1 | person, 
         ~ 1 + time | group, 
         emergence = ~ 1 + time, 
         method = "CEM2", 
         data = data6_miss)

summary.ce(m1)


