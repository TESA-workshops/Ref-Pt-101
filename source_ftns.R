#################################################################################################
# Survivorship
#################################################################################################

# default is f = 0, M can be a constant or vector of length n_ages
# n_ages is number of age classes from the age of recruitment to the plus group
# sel = vector for vulnerability-at-age and is of length (n_age) and only required if f > 0

survivorship_F <- function(f=0,M,n_ages,sel,message=T){
  if(n_ages<=2){message(paste0("number of age classes must be greater than 2"))}
  if(length(M)==1){
    if(message==T){message(paste0("Constant M = ",M," used for all ages"))}
    M<-rep(M,n_ages)
  }
  if(f==0){
    sel<-rep(0,n_ages)
    if(message==T){message("Assumed F = 0, unfished survivorship")}
  }
  if(n_ages != length(sel)){
    message("age clasees in vulnerability vector != n_ages")
    return(NA)
  }
  if(n_ages == length(sel)){
    l_age <- rep(NA,n_ages) 
    l_age[1] <- 1
    for(a in 2:(n_ages-1)){
      l_age[a] <- l_age[a-1]* exp(-(M[a-1]+f*sel[a-1]))
    }
    l_age[n_ages] <- l_age[n_ages-1]*exp(-(M[a-1]+f*sel[a-1]))/(1-exp(-(M[a]+f*sel[a])))
    return(l_age)
  }
}

#################################################################################################
# Reference Points
#################################################################################################

# M can be a constant or vector
# waa = weight-at-age vector
# mat = maturity-at-age vector
# sel = vulnerability-at-age vector
# all vectors must be of same length
# a and b are the Beverton_Holt a and b parameters
MSYcalc <- function(M,waa,mat,sel,a,b){
  
  output <- list()
  
  if(length(M)==1){
    message(paste0("Constant M = ",M," used for all ages"))
    M=rep(M,length(waa))
  }  
  if(length(waa) != length(mat) | length(waa) != length(sel)){
    message("at-age vectors are not the same length")
    return(NA)
  }
  
  f <- seq(0,5,0.001)
  DF <- data.frame(f=f, phi_f=NA, ypr_f=NA, eq_ssb_f=NA, eq_rec_f=NA, yield_f=NA)

    if(length(waa) == length(mat) & length(waa) == length(sel)){
    for(i in 1:length(f)){
      DF$phi_f[i] <- sum(survivorship_F(f=f[i],M=M,n_ages=length(sel),sel=sel,message=F)*waa*mat)
      DF$ypr_f[i] <- sum(survivorship_F(f=f[i],M=M,n_ages=length(sel),sel=sel,message=F)*waa*
                           (1-exp(-(M+f[i]*sel)))*f[i]*sel/(M+f[i]*sel))
    }  
    DF$eq_ssb_f <- (a*DF$phi_f-1)/b 
    DF$eq_ssb_f[DF$eq_ssb_f < 0] <- 0
    DF$eq_rec_f <- a*DF$eq_ssb_f/(1+b*DF$eq_ssb_f) 
    DF$eq_rec_f[DF$eq_rec_f < 0] <- 0
    DF$yield_f <- DF$ypr_f*DF$eq_rec_f
    output[[1]] <- f_msy <- f[which(DF$yield_f == max(DF$yield_f))] 
    output[[2]] <- msy <- DF$yield_f[which(DF$yield_f == max(DF$yield_f))]
    output[[3]] <- ssb_msy <- DF$eq_ssb_f[which(DF$yield_f == max(DF$yield_f))] 
    output[[4]] <- ssb_0 <- DF$eq_ssb_f[1]
    output[[5]] <- DF
    
    names(output) <- c("Fmsy","msy","SSBmsy","SSB0","DF")
    return(output)
  } 
}
