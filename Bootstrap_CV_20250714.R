# Example of CV evaluation under bootstrap, for one sample and one AIM, repeated data
# We compare the CV obtained after correction with  CV obtained after randomly selecting unstim values and applying same correction
# This gives the impact of the transformation, under no correlation between stim and unstim
# this is done within donor and compared to the CV of uncorrected stimulated data

# Gabrielle Boucher
# gabrielle.boucher@inflammgen.org
# Montreal Heart Institute
# July 14th 2025

CV <- function(V,...) return(sd(V,...)/mean(V,...))


x1 <- DAT_CMV.o.all[keep,k.aim]    # Stim for the donor "keep" and aim k.aim
x0 <- DAT_Unstim.o.all[keep,k.aim] # Unstim for the donor "keep" and aim k.aim

CV.raw <- CV(x1,na.rm=TRUE)
CV.0 <- CV( SI(x1, x0,0),na.rm=TRUE)
CV.1 <- CV( SI(x1, x0,1),na.rm=TRUE)
CV.05 <- CV( SI(x1, x0,0.5),na.rm=TRUE)
CV.05b <- CV( SI(x1, x0,0.5, corrected = FALSE),na.rm=TRUE)
CV.1b <- CV( SI(x1, x0,1, corrected = FALSE),na.rm=TRUE)


sd.boot.1 <- c()
sd.boot.0 <- c()
sd.boot.05 <- c()
sd.boot.05b <- c()
sd.boot.1b <- c()

mean.boot.1 <- c()
mean.boot.0 <- c()
mean.boot.05 <- c()
mean.boot.05b <- c()
mean.boot.1b <- c()
for(k in 1:100){
  s <- sample(1:9,9, replace=TRUE)
  sd.boot.1[k] <- sd( SI(x1, x0[s],1, corrected = TRUE),na.rm=TRUE)
  sd.boot.1b[k] <- sd( SI(x1, x0[s],1, corrected = FALSE),na.rm=TRUE)
  sd.boot.05[k] <- sd( SI(x1, x0[s],0.5, corrected = TRUE),na.rm=TRUE)
  sd.boot.05b[k] <- sd( SI(x1, x0[s],0.5, corrected = FALSE),na.rm=TRUE)
  sd.boot.0[k] <- sd( SI(x1, x0[s],0, corrected = TRUE),na.rm=TRUE)
  
  mean.boot.1[k] <- mean( SI(x1, x0[s],1, corrected = TRUE),na.rm=TRUE)
  mean.boot.1b[k] <- mean( SI(x1, x0[s],1, corrected = FALSE),na.rm=TRUE)
  mean.boot.05[k] <- mean( SI(x1, x0[s],0.5, corrected = TRUE),na.rm=TRUE)
  mean.boot.05b[k] <- mean( SI(x1, x0[s],0.5, corrected = FALSE),na.rm=TRUE)
  mean.boot.0[k] <- mean( SI(x1, x0[s],0, corrected = TRUE),na.rm=TRUE)
}

CV.boot.0 <- sd.boot.0/mean.boot.0
CV.boot.1 <- sd.boot.1/mean.boot.1
CV.boot.05 <- sd.boot.05/mean.boot.05
CV.boot.05b <- sd.boot.05b/mean.boot.05b
CV.boot.1b <- sd.boot.1b/mean.boot.1b


adjusted <- CV.raw/(c(CV.raw, median(CV.boot.0), median(CV.boot.1), median(CV.boot.05), median(CV.boot.05b),median(CV.boot.1b)) /
                      c(CV.raw, CV.0, CV.1, CV.05, CV.05b, CV.1b))


adjusted <-  c(CV.raw, CV.0, CV.1, CV.05, CV.05b, CV.1b) * 
  CV.raw/(c(CV.raw, median(CV.boot.0), median(CV.boot.1), median(CV.boot.05), median(CV.boot.05b),median(CV.boot.1b)) )
  


