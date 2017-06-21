v1 <- 3 #value
v2 <- .5
v3 <- 100

v1r <- c(1,5) #value # range
v2r <- c(0,1)
v3r <- c(75,125)

#column -> variables
#row -> scenarios

no.s <- 4 #no. of scenarios
vrl <- matrix(NA,3,3) #value range list
# vrl[1,] <- c(v1r,v1,no.s)
# vrl[2,] <- c(v2r, v2, no.s)
# vrl[3,] <- c(v3r, v3, no.s)
vrl[1,] <- c(v1,v1r)
vrl[2,] <- c(v2,v2r)
vrl[3,] <- c(v3,v3r)

simValueTable <- matrix(NA,nrow(vrl)*no.s,nrow(vrl)) #has row as scenario and column as variables
for(i in 1:nrow(vrl)){
  for(j in 1:(no.s)){
    simValueTable[j+(no.s*(i-1)),i] <- seq(vrl[i,2],vrl[i,3],length=no.s)[j]
  }
  simValueTable[is.na(simValueTable[,i]),i] <- vrl[i,1]
}

genSimValue <- function(vrl, no.s){
  simValueTable <- matrix(NA,nrow(vrl)*no.s,nrow(vrl)) #has row as scenario and column as variables
  for(i in 1:nrow(vrl)){
    for(j in 1:(no.s)){
      simValueTable[j+(no.s*(i-1)),i] <- seq(vrl[i,2],vrl[i,3],length=no.s)[j]
    }
    simValueTable[is.na(simValueTable[,i]),i] <- vrl[i,1]
  }
  simValueTable
}