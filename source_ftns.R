library(SAMtool) #Version 1.2.5
library(ggplot2)

getOMs<-function(n){
    In <- readRDS("In.rda")
    nsim <- In$OM@nsim
    proyears <- In$OM@proyears
    nyears <- In$OM@nyears
    nages <- In$OM@maxage+1
    hs <- c(0.8,0.6,0.95,rep(0.8,4))
    Ms <- c(rep(0.3,3),0.2,0.4,rep(0.3,2))
    mat <- In$OM@cpars$Mat_age[1,,nyears]
    v0 <- mat
    vl <- c(mat[2:12],1)
    vh <- c(0,mat[1:11])
    waa <- In$OM@cpars$Wt_age[1,,nyears]

        OMs <- list()
    
    for(i in 1:n){
      Inx <- In
      h <- hs[i]
      M <- Ms[i]
      if(i %in% 1:5) {sel <- v0}
      if(i == 6) {sel <- vl}
      if(i == 7) {sel <- vh}
      
      Inx$OM@h <- rep(h,2)
      Inx$OM@cpars$M_ageArray <- array(M,dim=c(nsim,nages,nyears+proyears))
      Inx$OM@cpars$V <- array(rep(sel,each=Inx$OM@nsim),dim=c(nsim,nages,nyears+proyears))
      Inx$OM@cpars$Mat_age <- array(rep(mat,each=Inx$OM@nsim),dim=c(nsim,nages,nyears+proyears))
      Inx$OM@cpars$Wt_age <- array(rep(NA,each=Inx$OM@nsim),dim=c(nsim,nages,nyears+proyears))
      Inx$OM@cpars$Wt_age[,,1:nyears] <- In$OM@cpars$Wt_age[,,1:nyears]
      Inx$OM@cpars$Wt_age[,,(nyears+1):(nyears+proyears)]  <- array(rep(In$OM@cpars$Wt_age[1,,nyears],each=Inx$OM@nsim),dim=c(nsim,nages,proyears))
      Inx$OM@cpars$Len_age <- array(rep(NA,each=Inx$OM@nsim),dim=c(nsim,nages,nyears+proyears))
      Inx$OM@cpars$Len_age[,,1:nyears] <- In$OM@cpars$Len_age[,,1:nyears]
      Inx$OM@cpars$Len_age[,,(nyears+1):(nyears+proyears)]  <- array(rep(In$OM@cpars$Len_age[1,,nyears],each=Inx$OM@nsim),dim=c(nsim,nages,proyears))
      
      Fit <- RCM_wrap(Inx)
      OMs[[i]]<-Fit
    }
    names(OMs)<-c("A","B","C","D","E","F","G")[1:n]
    return(OMs)
}

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
    M=rep(M,n_ages)
  }
  if(f==0){
    sel=rep(0,n_ages)
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
RPcalc <- function(M,waa,mat,sel,a,b){
  
  output <- list()
  
  if(length(M)==1){
    message(paste0("Constant M = ",M," used for all ages"))
    M=rep(M,length(waa))
  }  
  if(length(waa) != length(mat) | length(waa) != length(sel)){
    message("at-age vectors are not the same length")
    return(NA)
  }
  
  f <- seq(0,5,0.01)
  phi_f <- rep(NA,length(f))
  YPR_a <- SURV_a <- matrix(NA,ncol=length(waa),nrow=length(f))
  rownames(SURV_a)<-rownames(YPR_a)<-f
  
  if(length(waa) == length(mat) & length(waa) == length(sel)){
    for(i in 1:length(f)){
      SURV_a[i,] <- survivorship_F(f=f[i],M=M,n_ages=length(sel),sel=sel,message=F)
      for(j in 1:ncol(YPR_a)){
        YPR_a[i,j] <- SURV_a[i,j]*waa[j]*(1-exp(-(M[j]+f[i]*sel[j])))*f[i]*sel[j]/(M[j]+f[i]*sel[j])
      }
      phi_f[i] <- sum(SURV_a[i,]*waa*mat)
    }  
    ypr <- rowSums(YPR_a) 
    eq_rec_f <- (1/b*(a-1/phi_f))
    eq_rec_f[eq_rec_f<0] <- 0
    yield <- ypr*eq_rec_f
    output[[1]] <- f_msy <- f[which(yield==max(yield))]
    output[[2]] <- msy <- yield[which(yield==max(yield))]
    output[[3]] <- ssb_msy <- eq_rec_f[which(yield==max(yield))]/(a-eq_rec_f[which(yield==max(yield))]*b)
    output[[4]] <- ssb_0 <- max(eq_rec_f)/(a-max(eq_rec_f)*b)
    
    names(output) <- c("Fmsy","msy","SSBmsy","SSB0")
    return(output)
    
  } 
}


# Wrapper for the RCM fitting function that uses IN lists of data, OMs and misc settings

RCM_wrap<-function(IN){ # wrapper function allows alternative names for Input lists (IN)
  
  RCM(OM=IN$OM, data=IN$data, ESS=IN$misc$ESS, max_F=IN$misc$max_F,
      selectivity=IN$misc$selectivity, s_selectivity=IN$misc$s_selectivity,
      condition="catch2", mean_fit = T, rescale=1, LWT=IN$misc$LWT, cores=1,
      resample=T,comp_like=IN$misc$comp_like, map_log_rec_dev = c(1:(IN$OM@nyears-3), rep(NA, 3)),
      vul_par = IN$misc$vul_par, map_vul_par = IN$misc$map_vul_par, drop_nonconv = T)
  
}