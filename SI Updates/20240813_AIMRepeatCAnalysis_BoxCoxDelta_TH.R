library(openxlsx)
library(IDPmisc)
library(glue)
library(tidyverse)

original_data<-read.csv("AIM_Repeats_C_20230608.csv")
original_data$Repeat<- str_replace_all(original_data$Repeat, "-","_")


variables<-read.csv("Variables.csv")
variables<-variables %>% pull()

AIM_variables<-read.csv("AIM_variables.csv")
AIM_variables<-AIM_variables %>% pull()

Zero_to_NotZero<-function(x){
  ifelse((x<=0.01), 0.01, x)
}

original_data[7:37] <- lapply(original_data[7:37], Zero_to_NotZero)

SI <- function(x1,x0,L=0.5,H=NULL, s=1, e=0,OOB.V=1E-3, corrected=TRUE){
 
  # The function can be applied to vectors or matrix
  # x0 is unstim
  # x1 is stim
  
  # L is lambda
  # H is  Theta for correcting of the translation issue
  # You can set H as you wish, e.g set to 0 for small lambda and set to 1 for lambda close to 1
  # corrected =TRUE set a H value by default, which is F1(L) (see below)
 
  # setting H cancels the "corrected" option
  
  # OOB.V can be set to a low value, such as 1E-3
  # It can also be set to NaN or NA
  # This is the value returned when the SI is not defined; due to unstim >> stim
  
  # e and s are scaling factors. By default: no scaling
  
  
  F1 <- function(L){
    # S shaped function for the correction
    # returns 1 when L=1, 0 when L=0, 0.5 when L=0.5
    # with an S shape to get quick convergence to 0 and 1
    
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
  
  if(!is.null(H) & corrected){ corrected <- FALSE; warning("corrected set to FALSE\n")}
 
  if(corrected){
    H <- F1(L)
  }
  if(is.null(H)) H <- 0
  
  res <- NA*x1
  if(L>0 & L<=1){
    OOB <- which((x1^L+1-H)^(1/L) < x0 )   # out of bound (exclude NAs)
    N.OOB <- which((x1^L+1-H)^(1/L) >= x0 )   # not out of bound (exclude NAs)
    res[N.OOB] <- (x1[N.OOB]^L-x0[N.OOB]^L+1-H)^(1/L)
    res[OOB] <- OOB.V
  }
  if(L==0) res <- x1/x0
  
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

SI_data<-original_data %>% 
  select(Sample:CD3_EVENTS, CD4p_CD3p,Tregs_EVENTS,Tregs_CD4p,CD39p_EVENTS, CD39p_CD4p, CD8p_CD3p, MFI_CD4_CD4p, MFI_CD8_CD8p)

Stimulant <-(c("unstimulated", "CMV_peptides", "Infanrix", "CytoStim"))
DonorID <- (c("Z007-5", "Z029-4", "Z032-2", "Z048-2", "Z054-1"))

for (var in AIM_variables) {
  
  data<-original_data %>% 
    select(DonorID:DATE,{{var}}) %>% 
    mutate(across(DonorID:DATE, as.factor)) %>%
    mutate(y=(.data[[var]])) %>% 
    select(!all_of(var)) %>% 
    pivot_wider(names_from=Stimulant, values_from = y) %>% 
    mutate(bcxD_CMV_peptides=SI(x1 = CMV_peptides, x0 = unstimulated)) %>% 
    mutate(bcxD_Infanrix=SI(x1 = Infanrix, x0 = unstimulated)) %>% 
    mutate(bcxD_CytoStim=SI(x1 = CytoStim, x0 = unstimulated)) %>% 
    select(!(c(CMV_peptides, Infanrix, CytoStim))) %>% 
    pivot_longer(unstimulated:bcxD_CytoStim, names_to = "Stimulant", values_to = as.character(glue("{var}"))) %>% 
    mutate(across(6, round, 3))
  
  
  
  data$Stimulant<- str_replace_all(data$Stimulant, "bcxD_","")
  
  
  SI_data<-merge(SI_data, data, by=c("DonorID","Stimulant","Operator","Repeat","DATE"))
  print(var)
  
}
#just moving the Sample column to be the first column after failing to do so using other methods
SI_data<-SI_data %>% 
  relocate(Sample)

#saving all of the new subtracted values in a csv file.
write.csv(SI_data, "20240814_AIM Repeat C_cSIs.csv")

#to allow plotting on a log scale- add 0.001 to any 0 values
Zero_to_NotZero<-function(x){
  ifelse((x<=0.001), 0.001, x)
}




