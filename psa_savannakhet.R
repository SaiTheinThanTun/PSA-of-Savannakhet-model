v1 <- 3
v2 <- .5
v3 <- 100

v1r <- c(1,5)
v2r <- c(0,1)
v3r <- c(75,125)

#column -> variables
#row -> scenarios

no.s <- 4 #no. of scenarios
vrl <- matrix(NA,3,4)
vrl[1,] <- c(v1r,v1,no.s)
vrl[2,] <- c(v2r, v2, no.s)
vrl[3,] <- c(v3r, v3, no.s)

result <- matrix(NA,nrow(vrl)*no.s,nrow(vrl))
for(i in 1:nrow(vrl)){
  for(j in 1:(no.s)){
    result[j+(no.s*(i-1)),i] <- seq(vrl[i,1],vrl[i,2],length=no.s)[j]
  }
  result[is.na(result[,i]),i] <- vrl[i,3]
}