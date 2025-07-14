
############### CANASC Stimulation Index ###################
# Code to perform the box-cox transformation and compute the 
# Stimulation Index (SI)
# Written in base R only (no dependencies)

# Gabrielle Boucher 
# gabrielle.boucher@inflammgen.org
# Montreal Heart Institute
# July 14th 2025


# boxcox function (I added a default value of 0.5)
bc <- function(x,L=0.5){
  if(L==0) return(log(x)) else return((x^L-1)/L)
}
# inverse boxcox function
ibc <- function(x,L=0.5){
  if(L==0) return(exp(x)) else return((x*L+1)^(1/L))
}

# boxcox function modified for translation (I added a default value of 0.5)
bc.m <- function(x,L=0.5,H=NULL){
  F1 <- function(L){
    # S shaped function for the correction
    # returns 1 when L=1, 0 when L=0, 0.5 when L=0.5
    # with an Sigmoid shape to get quick convergence to 0 and 1
    
    sign0 <- function(v){
      s <- sign(v)
      s[s==0] <- 1
      return(s)
    } 
    
    L0 <- abs(L-0.5)+0.5
    f1 <- L0^((1-L0)/L0)
    
    f1 <- sign0(L-0.5)*f1 +(-sign0(L-0.5) +1)/2
    
    return(f1)
  }
  
  if(L==0) return(log(x)) 
  else{
    if(is.null(H)){
      H <- F1(L)
    }
    return((x^L-1)/L+H)
  } 
}

# to be noted; the box-cox is an estimation of log as lambda ->0, for values around 1
# unlike log, it can be computed for 0 (when lambda +0)
# most definition exclude x=0 because of the log case 

# for x=0 values, the box-cox will give: -1/lambda
# this is somehow equivalent to replacing 0s with exp(-1/lambda) in a log transformation
# eg at lambda = 0.5,  this is equivalent to replacing 0s with 0.14
# at lambda = 0.2, this is equivalent to replace 0s with 0.007

# One can do BOTH, ie you can add an epsilon before computing the box cox 
# With lambda < 0.1, the results rapidly converges to -Inf

# at lambda =1, boxcox is a linear transformation. We then have bc(x1,1) - bc(x0,1) = x1-x0
# If we get the inverse value from this, we get: x1-x0 +1
# the difference is then translated by 1

# Function SI below integrate the computing into a single function


SI <- function(x1,x0,L=0.5,H=NULL, s=1, e=0,OOB.V=1E-3, corrected=TRUE, rescale=FALSE){
  
  # The function can be applied to vectors or matrix
  # x0 is unstim
  # x1 is stim
  
  # L is lambda
  # H is  Theta for correcting of the translation issue
  # You can set H as you wish, e.g set to 0 for small lambda and set to 1 for lambda close to 1
  # corrected =TRUE set a H value by default, which is the sigmoid F1(L) (see below)
  
  # setting H cancels the "corrected" option
  
  # OOB.V can be set to a low value, such as 1E-3
  # It can also be set to NaN or NA
  # This is the value returned when the SI is not defined; due to unstim >> stim
  
  # e and s are scaling factors. By default: no scaling
  # e can be used to add a small epsilon to all values
  
  
  F1 <- function(L){
    # S shaped function for the correction
    # returns 1 when L=1, 0 when L=0, 0.5 when L=0.5
    # with an Sigmoid shape to get quick convergence to 0 and 1
    
    sign0 <- function(v){
      s <- sign(v)
      s[s==0] <- 1
      return(s)
    } 
    
    L0 <- abs(L-0.5)+0.5
    f1 <- L0^((1-L0)/L0)
    
    f1 <- sign0(L-0.5)*f1 +(-sign0(L-0.5) +1)/2
    
    return(f1)
  }
  
  
  d <- dim(x1)  # will be NULL if single value or vector, else will be nrow, ncol
  if(!is.null(d)) if(sum(dim(x0) != d)>0) stop("x0 and x1 not the same dimension\n")
  
  #scaling
  x1 <- as.vector(unlist((x1-e)/s))
  x0 <- as.vector(unlist((x0-e)/s))
  if(length(x0) != length(x1)) stop("x0 and x1 not the same length\n")
  
  if(!is.null(H) & corrected){ 
    corrected <- FALSE; 
    warning("corrected set to FALSE\n")
  }
  
  if(corrected){
    H <- F1(L)
  }
  if(is.null(H)) H <- 0
  
  res <- NA*x1
  if(L>0 & L<=1){
    OOB <- which((x1^L+1-L*H)^(1/L) < x0 )   # out of bound (exclude NAs)
    N.OOB <- which((x1^L+1-L*H)^(1/L) >= x0 )   # not out of bound (exclude NAs)
    res[N.OOB] <- (x1[N.OOB]^L-x0[N.OOB]^L+1-L*H)^(1/L)
    res[OOB] <- OOB.V
  }
  if(L==0) res <- x1/x0
  
  if(rescale){
    me1 <- mean(x1, na.rm=TRUE)
    sd1 <- sd(x1, na.rm=TRUE)
    res <- scale(res)*sd1+me1
  }
  
  dim(res) <- d
  return(res)
}


# Note: the sigmoid function F1 is also equivalent to:
F2 <- function(V){
  res <- NA*V
  w1 <- which(V<0.5)
  res[w1] <- 1-(1-V[w1])^( V[w1]/(1-V[w1]))
  w1 <- which(V>=0.5)
  res[w1] <- (V[w1])^( (1-V[w1])/V[w1])
  return(res)
}
# But the version F1 avoids some logical testing


# Apply SI to a vector of L values
A.SI <- function(x1, x0,LL, ...){
  RES <- matrix(NA,length(x1), length(LL))
  for(iL in 1:length(LL)){
    RES[,iL] <- SI(x1,x0,L=LL[iL], ...)
  }
  return(RES)
}
