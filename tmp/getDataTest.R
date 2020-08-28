## Test of getData() of Sherif data, 
## Date: 2020-08-28

library(MMP)
data(sherifdat)

# model on page 536 in Lang bookchapter (2019)
data <- getData(y ~ time, 
                ~ 1 | person, 
                ~ 1 | group, 
                emergence = ~ time, 
                method = "CEM", 
                data = sherifdat)

# model with interaction effect on page 536 in Lang bookchapter (2019) 
data <- getData(y ~ time, 
                ~ 1 | person, 
                ~ 1 | group, 
                emergence = ~ time * condition, 
                method = "CEM", 
                data = sherifdat)

### does not work, need to look at how to handle interaction effects