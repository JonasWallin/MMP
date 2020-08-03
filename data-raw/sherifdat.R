## code to prepare `DATASET` dataset goes here

## Creating Sherif data set

library(dplyr)
## Source: Muzafer Sherif, "A study of some social factors in perception: Chapter 3."
## Archives of Psychology, 1935, 27, No. 187,  23- 46

sherifdat=expand.grid(person = 1:3,time = -1:3,
                      group = 1:8)

sherifdat$y <- c(## Starting with individual session
                  # group 1
                    7.50, 0.66, 1.83,
                    3.90, 2.01, 2.75,
                    3.78, 2.11, 2.64,
                    2.10, 2.18, 2.29,
                    NA, NA, NA,
                  # group 2
                    1.92, 0.14, 1.42,
                    1.31, 1.03, 1.49,
                    1.94, 1.45, 1.68,
                    1.79, 1.45, 1.58,
                    NA, NA, NA,
                  # group 3
                    1.84, 2.66, 4.26,
                    1.76, 1.85, 2.63,
                    2.17, 2.25, 2.83,
                    1.79, 1.43, 1.79,
                    NA, NA, NA,
                  # group 4
                    2.20, 4.94, 3.57,
                    3.51, 3.55, 3.36,
                    4.69, 4.56, 3.79,
                    4.29, 4.62, 4.35,
                    NA, NA, NA,
                ## Starting with the group situation
                  # group 5
                    NA, NA, NA,
                    4.30, 3.65, 4.29,
                    4.13, 3.91, 4.03,
                    3.92, 3.52, 3.77,
                    4.39, 3.41, 3.58,
                  # group 6
                    NA, NA, NA,
                    4.55, 3.96, 3.88,
                    1.85, 1.58, 1.80,
                    1.39, 1.63, 1.46,
                    2.00, 1.18, 1.39,
                  # group 7
                    NA, NA, NA,
                    2.55, 2.52, 2.22,
                    1.60, 2.26, 1.76,
                    1.79, 1.95, 1.66,
                    1.85, 1.58, 2.38,
                  # group 8
                    NA, NA, NA,
                    1.44, 2.54, 3.17,
                    5.25, 5.23, 5.83,
                    4.31, 4.50, 4.47,
                    4.14, 4.70, 4.45
                  )

# condition = 1 if Starting with individual session, = 0 if Starting with the group situation
sherifdat$condition<-rep(c(1,0),each=60)


## Centering
# The within-group means for each time point are
# subtracted from the observed y values.

means <- sherifdat %>%
  group_by(group, time) %>%
  summarise_at(vars(y), list(~mean(., na.rm=TRUE)))

sherifdat <- left_join(sherifdat, means, by =c("group", "time"))

sherifdat$y.centered <- sherifdat$y.x-sherifdat$y.y

sherifdat <- ungroup(sherifdat)

sherifdat <- rename(sherifdat, y = y.x, g.mean = y.y)

sherifdat <- na.omit(sherifdat)

usethis::use_data(sherifdat, sherifdat, overwrite = TRUE)
