source("charma_fit.r")
## ============================================================================================

# function for generating Chen observations using the inversion method   
rchen<-function(n,mu,lambda){
  u=runif(n)
  y= (log(1-(log(1-u)/(log(1-0.5)/(1-exp(mu^lambda)))))) ^(1/lambda)
  return(y)
}

simu.CHARMA <- function(n,beta0,lambda,phi=NA,theta=NA,beta=NA,freq=12,link="log") {
  ar <- NA
  ma <- NA
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "identity")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"log\"  and \"identity\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  ###### ARMA model
  if(any(is.na(phi)==F) && any(is.na(theta)==F))
  {
    p <- max(ar)
    q <- max(ma)
    m <- 50
    maxx=max(p,q)
    
    ynew <-rep(beta0,(n+m))
    mu <- linkinv(ynew)
    
    error <- rep(0,n+m) # E(error)=0 
    eta <- y <- NULL
    
    for(i in (maxx+1):(n+m))
    {
      eta[i] <- beta0 + (phi%*%(ynew[i-ar])) + (theta%*%error[i-ma]) 
      mu[i] <- linkinv(eta[i])
      y[i] <- rchen(1,mu=mu[i],lambda) 
      ynew[i] <- linkfun(y[i]) 
      error[i] <- ynew[i]-eta[i]   
    }
    return(ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  # AR model
  if(any(is.na(phi)==F) && any(is.na(theta)==T))
  {
    p <- max(ar)
    m <- 50
    
    ynew <- rep(beta0,(n+m))
    mu <- linkinv(ynew)
    
    eta <- y <- NULL
    
    eta <- beta0 + (phi%*%(ynew[(p+1)-ar]))
    muu <- linkinv(eta)
    
    for(i in (p+1):(n+m))
    {
      eta[i] <- beta0 + (phi%*%(ynew[i-ar]))
      mu[i] <- linkinv(eta[i])
      y[i] <- rchen(1,mu=mu[i],lambda)
      ynew[i] <- linkfun(y[i]) 
      
    }
    return(ts(y[(m+1):(n+m)],frequency=freq) )
  } 
  
  # MA model
  if(any(is.na(phi)==T) && any(is.na(theta)==F))
  {
    q <- max(ma)
    m <- 50
    
    ynew <-rep(beta0,(n+m))
    mu <- linkinv(ynew)
    
    eta <- y <- error <- rep(0,n+m) 
    X <- cbind(sin(2*pi*(1:(n+m))/50))
    
    for(i in (q+1):(n+m))
    {
      eta[i] <- beta0 + (theta%*%error[i-ma]) 
      mu[i] <- linkinv(eta[i])
      y[i] <- rchen(1,mu=mu[i],lambda) 
      ynew[i] <- linkfun(y[i])  
      error[i]<- ynew[i]-eta[i]
    }
    return(ts(y[(m+1):(n+m)],frequency=freq) )
  } 
}

## ===========================================================================================
# Monte Carlo simulation - example for CHARMA(1,0)
vn <-c(100, 250, 500)
lambda <-c(0.7)
beta0 <- c(0.2)
phi <- c(0.3)

R <- 1000
k <- 3

for(n in vn)
{
  cont=0
  mreplicas=matrix(rep(NA, R*k), ncol=k)
  i=1
  vetory=NULL
  while(i<=R)
  {
    data_simu <- simu.CHARMA(n,beta0,lambda,phi,freq=12,link="log")
    fit <- try(charma.fit(y=data_simu,ar=c(1),h1=0, graf=0,desc=0), T)
    
    if(class(fit) == "try-error"){cont=cont+1} 
    if(class(fit) != "try-error")
    {
      
      if(fit$conv!=0){cont=cont+1}
      mreplicas[i,]=c(fit$coeff)
      
      if((100*i/R)%%10 == 0)
        print(c(100*i/R, "%"), quote = F)
      i=i+1
      
    }
  }
  
  par <- c(beta0,phi,lambda)
  mean <- colMeans(mreplicas)
  bias <- mean - c(par)
  SD <- apply(mreplicas,2, sd)
  MSE <- (bias^2)+(SD^2)
  RB <- (bias/c(par))*100
  
  mresults <- rbind(mean, bias, RB, SD, MSE)
  
  colnames(mresults) <-c("beta0","phi","lambda")
  rownames(mresults) <-c("mean","bias","RB(%)","SD","MSE")
  print(" ", quote=F)
  print(c("Sample sizes =",n),quote=F) 
  print(round(mresults,3))
  
}

## ===========================================================================================
# Monte Carlo simulation - example for CHARMA(0,1)
vn <-c(100, 250, 500)
lambda <-c(0.7)
beta0 <- c(0.2)
theta <- c(0.2)

R <- 1000
k <- 3

for(n in vn)
{
  cont=0
  mreplicas=matrix(rep(NA, R*k), ncol=k)
  i=1
  vetory=NULL
  while(i<=R)
  {
    data_simu <- simu.CHARMA(n,beta0,lambda,theta,freq=12,link="log")
    fit <- try(charma.fit(y=data_simu,ma=c(1),h1=0, graf=0,desc=0), T)
    
    if(class(fit) == "try-error"){cont=cont+1} 
    if(class(fit) != "try-error")
    {
      
      if(fit$conv!=0){cont=cont+1}
      mreplicas[i,]=c(fit$coeff)
      
      if((100*i/R)%%10 == 0)
        print(c(100*i/R, "%"), quote = F)
      i=i+1
      
    }
  }
  
  par <- c(beta0,theta,lambda)
  mean <- colMeans(mreplicas)
  bias <- mean - c(par)
  SD <- apply(mreplicas,2, sd)
  MSE <- (bias^2)+(SD^2)
  RB <- (bias/c(par))*100
  
  mresults <- rbind(mean, bias, RB, SD, MSE)
  
  colnames(mresults) <-c("beta0","theta","lambda")
  rownames(mresults) <-c("mean","bias","RB(%)","SD","MSE")
  print(" ", quote=F)
  print(c("Sample sizes =",n),quote=F) 
  print(round(mresults,3))
  
}


## ===========================================================================================
# Monte Carlo simulation - example for CHARMA(1,1)
vn <-c(100, 250, 500)
lambda <-c(0.7)
beta0 <- c(0.3)
phi <- c(0.2)
theta <- c(0.3)

R <- 1000
k <- 4

for(n in vn)
{
  cont=0
  mreplicas=matrix(rep(NA, R*k), ncol=k)
  i=1
  vetory=NULL
  while(i<=R)
  {
    data_simu <- simu.CHARMA(n,beta0,lambda,phi,theta,freq=12,link="log")
    fit <- try(charma.fit(y=data_simu,ar=c(1),ma=c(1), h1=0, graf=0,desc=0), T)
    
    if(class(fit) == "try-error"){cont=cont+1} 
    if(class(fit) != "try-error")
    {
      
      if(fit$conv!=0){cont=cont+1}
      mreplicas[i,]=c(fit$coeff)
      
      if((100*i/R)%%10 == 0)
        print(c(100*i/R, "%"), quote = F)
      i=i+1
      
    }
  }
  
  par <- c(beta0,phi,theta,lambda)
  mean <- colMeans(mreplicas)
  bias <- mean - c(par)
  SD <- apply(mreplicas,2, sd)
  MSE <- (bias^2)+(SD^2)
  RB <- (bias/c(par))*100
  
  mresults <- rbind(mean, bias, RB, SD, MSE)
  
  colnames(mresults) <-c("beta0","phi","theta","lambda")
  rownames(mresults) <-c("mean","bias","RB(%)","SD","MSE")
  print(" ", quote=F)
  print(c("Sample sizes =",n),quote=F) 
  print(round(mresults,3))
  
}