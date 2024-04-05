#Load libraries
#Note: tidyverse contains ggplot2, dplyr, tidyr, stringr
library(openxlsx)
library(IDPmisc)
library(glue)
library(tidyverse)

#this is just to create a set of graphs without any parameters highlighted
all_data<-read.csv("20240316_RawPREVENTDataForBoxcoxAnalysis.csv")
variables<-read.csv("AIM_variables.csv")
variables<-variables %>% pull()
stimulants<-c("DMSO","Fluzone", "COVID_WT", "COVID_BA4_5")

Neg_to_Zero<-function(x){
  ifelse((x<=0.005), 0.005, x)
}

# boxcox function (I added a default lambda value of 0.5)
#default lambda value of 0.5
bc <- function(x,l=0.5){if(l==0) return(log(x)) else return((x^l-1)/l)}
# inverse boxcox function
ibc <- function(x,l=0.5){if(l==0) return(exp(x)) else return((x*l+1)^(1/l))}

#var<-"CD4_CD134pCD25p"

all_data <- all_data %>%
  filter(Stim != "Fluzone" & Stim != "Cytostim") %>%
  filter(DonorID != "UBC_004" & DonorID != "UBC_059") %>%
  filter(Timepoint != "V1" & Timepoint != "6M" & Timepoint != "1Y")


all_data <- all_data %>%
  arrange(Stim, Timepoint, DonorID)
all_data

for (var in variables) {
  
  subtracted_data<-all_data %>% 
    select(DonorID:Stim,{{var}}) %>% 
    pivot_wider(names_from = Stim, values_from = {{var}}) %>% 
    mutate(Covid_WT_sum=COVID_WT-DMSO, Covid_BA4_5_sum=COVID_BA4_5-DMSO) %>% 
    pivot_longer(DMSO:Covid_BA4_5_sum, names_to="Stimulant", values_to = {{var}}) %>%
    mutate(across(DonorID:Stimulant,as.factor)) %>%
    filter(Stimulant != "COVID_WT" & Stimulant != "COVID_BA4_5") %>%
    mutate(Stimulant=fct_recode(Stimulant, COVID_WT="Covid_WT_sum", COVID_BA4_5="Covid_BA4_5_sum")) %>% 
    mutate(Timepoint=fct_relevel(Timepoint,"VY","V2","V3"))
  
  #use lapply to apply the Neg to Zero function to all rows
  #subtracted_data[[var]] <- lapply(subtracted_data[[var]], Neg_to_Zero)
  #subtracted_data[[var]] <- unlist(subtracted_data[[var]]) 
  subtracted_data <- arrange(subtracted_data, Stimulant, Timepoint, DonorID)
  write.csv(subtracted_data, glue("PREVENT_Subtracted_{var}.csv"))
  
}

for (var in variables) {
  
  bcxsubtracted_data<-all_data %>% 
    select(DonorID:Stim,{{var}}) %>% 
    pivot_wider(names_from = Stim, values_from = {{var}}) %>% 
    mutate(Covid_WT_sum=ibc(bc(COVID_WT)-bc(DMSO)), Covid_BA4_5_sum=ibc(bc(COVID_BA4_5)-bc(DMSO))) %>% 
    pivot_longer(DMSO:Covid_BA4_5_sum, names_to="Stimulant", values_to = {{var}}) %>%
    mutate(across(DonorID:Stimulant,as.factor)) %>%
    filter(Stimulant != "COVID_WT" & Stimulant != "COVID_BA4_5") %>%
    mutate(Stimulant=fct_recode(Stimulant, COVID_WT="Covid_WT_sum", COVID_BA4_5="Covid_BA4_5_sum")) %>% 
    mutate(Timepoint=fct_relevel(Timepoint,"VY","V2","V3"))
  
  #use lapply to apply the Neg to Zero function to all rows
  #bcxsubtracted_data[[var]] <- lapply(bcxsubtracted_data[[var]], Neg_to_Zero)
  #bcxsubtracted_data[[var]] <- unlist(bcxsubtracted_data[[var]]) 
  bcxsubtracted_data <- arrange(bcxsubtracted_data, Stimulant, Timepoint, DonorID)
  write.csv(bcxsubtracted_data, glue("PREVENT_Boxcox_{var}.csv"))
}

