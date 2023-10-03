# How they calculate ICC1 in multilevel package (by Paul Bliese)

library(multilevel)

# use formula (1) from Bartko (1976), On various intraclass correlation reliability coefficients,
# Psychological Bulletin Vol 83, No 5, p. 762-765

data(bh1996)
# one way anova using lm function:
hrs.mod<-aov(HRS~as.factor(GRP),data=bh1996)
ICC1(hrs.mod)

# what it does:
MOD <- summary(hrs.mod)
MSB <- MOD[[1]][1, 3] # mean sum of squares between groups
MSW <- MOD[[1]][2, 3] # mean residual sum of squares
GSIZE <- (MOD[[1]][2, 1] + (MOD[[1]][1, 1] + 1))/(MOD[[1]][1,1] + 1) # group size
# assumes a balanced design (=the same number of observations in each group)
OUT <- (MSB - MSW)/(MSB + ((GSIZE - 1) * MSW))


