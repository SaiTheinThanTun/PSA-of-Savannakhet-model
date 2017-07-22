#permutation testing
library(gtools)

permutations(3,2,c('l','m','h'),repeats.allowed = TRUE)
permutations(3,2,1:3,repeats.allowed = TRUE)
