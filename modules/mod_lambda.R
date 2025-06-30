library(shiny)
library(MASS)
library(tidyverse)


# Functions from given code
bc <- function(x,L=0.5){if(L==0) return(log(x)) else return((x^L-1)/L)}
ibc <- function(x,L=0.5){if(L==0) return(exp(x)) else return((x*L+1)^(1/L))}

test.lambda.lik <- function(x1, x0, LL=NULL,f1=NULL, method="rlm", DAT=NULL, plotit=TRUE){
  # "best" lambda here will have slope = 0 on this model f(x1)-f(x0) ~ f(x0) 
  # (equivalent to slope of 1 on f(x1)~f(x0) )
  # We use the density distribution of the slope estimate as the posterior distribution
  # under a non informative prior
  # the distribution is non-central t we can approximate by normal
  # normal would put more weight on the estimated coef
  # regression is sensitive to extreme values and distribution should be considered an approximation
  # rlm is less sensitive to outliers
  # we assume independence of the samples
  # Spearman rho's is a more robust method, but we do not compute a posterior
  # here we simply return the rho value, which should be ~0 at best lambda
  
  # Note: negatives estimate for rho correspond to "over-corrected"
  # if x1 and x0 are not related, we expect negatives rho
  # if x1 and x0 are unrelated, the maximum for post.L is expected to be close to 1
  # when there is no lambda with beta=0, there is no reliable estimate of lambda 
  
  
  # This version also incorporates the likelihood with a fixed beta
  # data is corrected for the Jacobian of the transformation
  
  # f1 accounts for covariates. This changes the interpretation of the results
  
  if(is.null(LL)) LL <- seq(0,1,0.01)
  nona <- !is.na(x0)&!is.na(x1)
  x0 <- x0[nona]
  x1 <- x1[nona]
  n <- length(nona)
  
  DAT <- DAT[nona,]
  
  if(is.null(f1)){
    fb1 <- as.formula("yy~x")
    fb2 <- as.formula("y~offset(x)")
  }
  if(!is.null(f1)){
    #consider covariate from f1
    fb1 <- as.formula(paste0("yy~x+", f1) )
    fb2 <- as.formula(paste0("y~offset(x)+", f1) )
  }
  
  post.b <- c()
  rho <- c()
  beta <- c()
  L.Lik <- c()
  
  for(i in 1:length(LL)){
    JJ <- (LL[i]-1)*sum(log(x1)) # log Jacobian
    K <- (exp(JJ))^(1/n) #data rescaling factor to account for Jacobian
    
    y <- bc(x1, LL[i])
    x <- bc(x0, LL[i])
    
    #We rescale both x and y to keep the relationship
    y <- y/K
    x <- x/K
    yy <- y-x
    
    if(method=="lm"){
      Mb <- lm(fb1,data=DAT)
      M1 <- lm(fb2,data=DAT)
    }
    if(method=="rlm"){
      Mb <- rlm(fb1,data=DAT)
      M1 <- rlm(fb2,data=DAT)
    }
    beta[i] <- Mb$coef[2]
    resid1 <- resid(M1)
    rho[i] <- cor( resid1,x, method="spearman", use="complete.obs")
    #dfs <- Mb$df.residual
    dfs <- n - length(Mb$coefficients) # for rlm compatibility
    post.b[i] <- dt((0-Mb$coef[2])/summary(Mb)$coef[2,2],dfs)
    #post.b[i] <- dt((1-Mb$coef[2])/summary(Mb)$coef[2,2],dfs)
    
    Sig2 <- sum(resid1^2)/n
    tmp <- anova(M1)
    #J <- (LL[i]-1)*sum(log(x1)) # log Jacobian if not corrected data
    J <- 0 # if jacobian already included in data transformation
    L.Lik[i] <- (-n/2)*log(Sig2)+J
    #L.Lik[i] <- logLik(M1)+J   # equivalent
    
  }
  
  # under non informative prior for L, beside being 0-1 (here being LL)
  # we get the posterior P(L | x1,x0,b=0)
  # post.L <- post.b/sum(post.b)
  post.L <- post.b/sum(post.b, na.rm=TRUE)
  Lik <- exp(L.Lik-median(L.Lik))
  Lik <- Lik/sum(Lik)
  
  res <- data.frame(LL=LL,post.L=post.L, beta, rho, L.Lik, Lik)
  
  if(plotit){
    tmp<- plotres(res)
    print(tmp)
  }
  
  return(res)
}


plotres <- function(MyRes){
  tolerance <- 0.1
  L.est <- data.frame(method = c("beta", "post.L", "rho", "L.Lik", "L.Lik.LB", "L.Lik.UB"), L=NA)
  
  with(MyRes, {
    w0 <- which.min(abs(beta))
    par(mfrow=c(2,2))
    plot(LL,beta, type="l")
    abline(h=0)
    if( abs(beta[w0])< tolerance ){
      abline(v=LL[w0], col="blue")
      text(LL[w0], beta[w0],paste("L=",round(LL[w0],2)), pos=2, col="blue" )
      L.est$L[1] <<- LL[w0]
    }
    
    w1 <- which.max(post.L)
    plot(LL,post.L, type="l")
    if( abs(beta[w1])< tolerance ){
      abline(v=LL[w1], col="blue")
      text(LL[w1], post.L[w1],paste("L=",round(LL[w1],2)), pos=2, col="blue" )
      L.est$L[2] <<- LL[w1]
    }
    plot(LL, rho, type="l", ylim=c(-0.5, 0.5))
    abline(h=0, col="red")
    w2 <- which.min(abs(rho))
    if( abs(rho[w2])< tolerance ){
      abline(v=LL[w2],col="red")
      text(LL[w1], rho[w2],paste("L=",round(LL[w2],2)), pos=2, col="red" )
      L.est$L[3] <<- LL[w2]
    }
    w3 <- which.max(L.Lik)
    plot(LL,L.Lik, type="l")
    abline(v=LL[w3], col="purple")
    text(LL[w3], L.Lik[w3],paste("L=",round(LL[w3],2)), pos=2, col="purple" )
    L.est$L[4] <<- LL[w3]
    
    mL <- L.Lik[w3]
    delta <- mL-1.92
    w.ub.l <- which(LL > LL[w3])
    w.ub <- w.ub.l[which.min( abs(L.Lik-delta)[w.ub.l])]
    w.lb.l <- which(LL < LL[w3])
    w.lb <- w.lb.l[which.min( abs(L.Lik-delta)[w.lb.l])]
    #abline(h=mL, col="purple")
    abline(h=delta, col="purple", lty=2)
    if(length(w.ub)>0){
      if( abs(L.Lik[w.ub]-delta)< tolerance*abs(mL) ){
        abline(v=LL[w.ub], col="purple", lty=2)
        L.est$L[6] <<- LL[w.ub]
      }}
    if(length(w.lb)>0){
      if( abs(L.Lik[w.lb]-delta)< tolerance*abs(mL) ){
        abline(v=LL[w.lb], col="purple", lty=2)
        L.est$L[5] <<- LL[w.lb]
      }}
  })
  return(L.est)
}

plotres2 <- function(MyRes){
  tolerance <- 0.1
  L.est <- data.frame(method = c("beta", "post.L", "rho", "L.Lik", "L.Lik.LB", "L.Lik.UB"), L=NA)
  
  with(MyRes, {
    
    w3 <- which.max(L.Lik)
    plot(LL,L.Lik, type="l")
    abline(v=LL[w3], col="purple")
    text(LL[w3], L.Lik[w3],paste("L=",round(LL[w3],2)), pos=2, col="purple" )
    L.est$L[4] <<- LL[w3]
    
    mL <- L.Lik[w3]
    delta <- mL-1.92
    w.ub.l <- which(LL > LL[w3])
    w.ub <- w.ub.l[which.min( abs(L.Lik-delta)[w.ub.l])]
    w.lb.l <- which(LL < LL[w3])
    w.lb <- w.lb.l[which.min( abs(L.Lik-delta)[w.lb.l])]
    #abline(h=mL, col="purple")
    abline(h=delta, col="purple", lty=2)
    if(length(w.ub)>0){
      if( abs(L.Lik[w.ub]-delta)< tolerance*abs(mL) ){
        abline(v=LL[w.ub], col="purple", lty=2)
        L.est$L[6] <<- LL[w.ub]
      }}
    if(length(w.lb)>0){
      if( abs(L.Lik[w.lb]-delta)< tolerance*abs(mL) ){
        abline(v=LL[w.lb], col="purple", lty=2)
        L.est$L[5] <<- LL[w.lb]
      }}
  })
  return(L.est)
}


plotres3 <- function(MyRes){
  tolerance <- 0.1
  L.est <- data.frame(method = c("beta", "post.L", "rho", "L.Lik", "L.Lik.LB", "L.Lik.UB"), L=NA)
  
  with(MyRes, {
    w0 <- which.min(abs(beta))
    #par(mfrow=c(2,2))
    plot(LL,beta, type="l")
    abline(h=0)
    if( abs(beta[w0])< tolerance ){
      abline(v=LL[w0], col="blue")
      text(LL[w0], beta[w0],paste("L=",round(LL[w0],2)), pos=2, col="blue" )
      L.est$L[1] <<- LL[w0]
    }
    
    w1 <- which.max(post.L)
    plot(LL,post.L, type="l")
    if( abs(beta[w1])< tolerance ){
      abline(v=LL[w1], col="blue")
      text(LL[w1], post.L[w1],paste("L=",round(LL[w1],2)), pos=2, col="blue" )
      L.est$L[2] <<- LL[w1]
    }
    plot(LL, rho, type="l", ylim=c(-0.5, 0.5))
    abline(h=0, col="red")
    w2 <- which.min(abs(rho))
    if( abs(rho[w2])< tolerance ){
      abline(v=LL[w2],col="red")
      text(LL[w1], rho[w2],paste("L=",round(LL[w2],2)), pos=2, col="red" )
      L.est$L[3] <<- LL[w2]
    }
    w3 <- which.max(L.Lik)
    plot(LL,L.Lik, type="l")
    abline(v=LL[w3], col="purple")
    text(LL[w3], L.Lik[w3],paste("L=",round(LL[w3],2)), pos=2, col="purple" )
    L.est$L[4] <<- LL[w3]
    
    mL <- L.Lik[w3]
    delta <- mL-1.92
    w.ub.l <- which(LL > LL[w3])
    w.ub <- w.ub.l[which.min( abs(L.Lik-delta)[w.ub.l])]
    w.lb.l <- which(LL < LL[w3])
    w.lb <- w.lb.l[which.min( abs(L.Lik-delta)[w.lb.l])]
    #abline(h=mL, col="purple")
    abline(h=delta, col="purple", lty=2)
    if(length(w.ub)>0){
      if( abs(L.Lik[w.ub]-delta)< tolerance*abs(mL) ){
        abline(v=LL[w.ub], col="purple", lty=2)
        L.est$L[6] <<- LL[w.ub]
      }}
    if(length(w.lb)>0){
      if( abs(L.Lik[w.lb]-delta)< tolerance*abs(mL) ){
        abline(v=LL[w.lb], col="purple", lty=2)
        L.est$L[5] <<- LL[w.lb]
      }}
  })
  return(L.est)
}

# UI and Server
lambdaUI <- function(id, title = "Lambda (λ) Estimation Tools") {
  ns <- NS(id)

  tabPanel(
    title,
    sidebarPanel(
      fileInput(ns("csv"), "Upload CSV", accept = ".csv"),
      textInput(ns("stim_col"),   "Stimulant Column Name",  "Stim"),
      textInput(ns("stim_lvl"),   "Stimulant",            "COVID_WT"),
      textInput(ns("unstim_lvl"), "Unstimulated Parameter",          "DMSO"),
      numericInput(ns("eps"), "Epsilon", value = 0.001, step = 0.0005),
      uiOutput(ns("aim_selector")),
      actionButton(ns("run"), "Run")
    ),
    mainPanel(
      plotOutput(ns("lambda_plot"), height = 600)
    )
  )
}


lambdaServer <- function(id) {
  moduleServer(id, function(input, output, session) {

    # populate AIM variable selector once the CSV is in
    output$aim_selector <- renderUI({
      req(input$csv)
      df <- read.csv(input$csv$datapath, nrows = 1, stringsAsFactors = FALSE)
      selectInput(session$ns("aim_var"), "AIM Variable (column)",
                  choices = names(df), selected = names(df)[1])
    })

    # run on "Run"
    observeEvent(input$run, {
      req(input$csv, input$aim_var, input$stim_col,
          nzchar(input$stim_lvl), nzchar(input$unstim_lvl))

      # User input
      epsilon <- input$eps
      data    <- read.csv(input$csv$datapath, stringsAsFactors = FALSE)

      Stim   <- data %>% filter(.data[[input$stim_col]] == input$stim_lvl)
      Unstim <- data %>% filter(.data[[input$stim_col]] == input$unstim_lvl)

      Stim   <- Stim[[input$aim_var]]   + epsilon
      Unstim <- Unstim[[input$aim_var]] + epsilon

      # open a fresh graphics device inside renderPlot
      output$lambda_plot <- renderPlot({
        test.lambda.lik(x1 = Stim, x0 = Unstim)  # this ultimately prints the plots
      })
    })
  })
}
