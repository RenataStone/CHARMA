
charma.fit<- function (y,ar=NA,ma=NA,link="log",h1=0,X=NA,X_hat=NA,graf=0,desc=1)
{
  maxit1<-100
  tau <- 0.5

  z <- c()
  link <- make.link("log")
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- function(t) 1/(link$mu.eta(link$linkfun(t)))
  
  y <- as.vector(y)
  ynew <- linkfun(y)
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  dchen <- function (x, b=1, lambda=1, log=FALSE)
  {
    pdf=x
    pdf[log==FALSE]=b*lambda*x**(b-1)*exp(x**b)*exp(lambda-lambda*exp(x**b))
    pdf[log==TRUE]=log(b*lambda)+(b-1)*log(x)+x**b+lambda-lambda*exp(x**b)
    return(pdf)
  }
  
  pchen <- function (x, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
  {
    cdf=x
    cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(lambda-lambda*exp(x**b))
    cdf[log.p==FALSE&lower.tail==FALSE]=exp(lambda-lambda*exp(x**b))
    cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(lambda-lambda*exp(x**b)))
    cdf[log.p==TRUE&lower.tail==FALSE]=lambda-lambda*exp(x**b)
    return(cdf)
  }
  
  y_prev <- c(rep(NA,(n+h1)))

  # initial values
  if(any(is.na(ar)==F)) # with AR
  {
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar]
    }
    
    Z <- cbind(rep(1,(n-m)),P)
  }else{
    Z <- as.matrix(rep(1,(n-m)))
  }
  
  if(any(is.na(X)==T)) # without regressors
  {
    x <- as.matrix(Z)
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    lambda<- 0.7 ### lambda
    
  }else{
    X_hat <- as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    #Ystar = log(Y/(1-Y))
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((n1 - k) * (dlink)^2)  ###
    lambda<- 0.7 ### lambda
  }
  ############
  ######### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T))
  {
    reg <- c(mqo, rep(0,q1), lambda) # initializing the parameter values
    
    loglik <- function(z) 
    {
      beta0 <- z[1]
      phi <- z[2:(p1+1)] 
      theta <- z[(p1+2):(p1+q1+1)]
      lambda <- z[p1+q1+2] # precision parameter
      
      error <- rep(0,n) # E(error)=0 
      eta <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
      
      ETA <- (log(1- tau))/(1-(exp(mu^lambda))) 
      ll <- suppressWarnings( log( dchen(y1, b=lambda, lambda=ETA) ) )
      sum(ll)
    } 
    
    escore <- function(z)
    {
      beta0 <- z[1]
      phi <- z[2:(p1+1)] 
      theta <- z[(p1+2):(p1+q1+1)]
      lambda <- z[p1+q1+2] # precision parameter
      
      error <- rep(0,n) # E(error)=0 
      eta <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n] 
      
      ######
      Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        Q[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrences
      deta.dbeta0 <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta <- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
        deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }
      
      b0 <- deta.dbeta0[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rQ <- deta.dtheta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      v <- ((lambda*exp(mu^lambda)*mu^(lambda-1)))/(1-exp(mu^lambda))+
        ((lambda*log(0.5)*exp(mu^lambda)*(mu^(lambda-1))*(1-exp(y1^lambda)))/((1-exp(mu^lambda))^2))
      
      c <- ((mu^(lambda)*exp(mu^lambda)*log(mu))/(1-exp(mu^lambda))) + (1/lambda)+
        ((log(0.5)*exp(mu^lambda)*(mu^lambda)*log(mu)*(1-exp(y1^lambda)))/(1-exp(mu^lambda))^2)+
        log(y1) - ((log(0.5)*exp(y1^lambda)*(y1^lambda)*log(y1))/(1-exp(mu^lambda))) + ((y1^lambda)*log(y1))
      
      Ubeta0 <- t(b0) %*% mT %*% v
      Uphi <-   t(rP) %*% mT %*% v
      Utheta <- t(rQ) %*% mT %*% v
      Ulambda <- sum(c) #t(c)              
      
      rval <- c(Ubeta0,Uphi,Utheta,Ulambda)
      return(rval)
      
    }
    names_phi <- c(paste("phi",ar,sep=""))
    names_theta <- c(paste("theta",ma,sep=""))
    names_par <- c("beta0",names_phi,names_theta,"lambda") 
    
    opt <- optim(reg, loglik, gr=escore, hessian=T, method = "BFGS", 
                 control = list(fnscale = -1))
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+2)]
    names(coef) <- names_par
    z$coeff <- coef
    
    beta0 <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    lambda <- coef[p1+q1+2] # precision parameter
    
    z$beta0 <- beta0
    z$phi <- phi
    z$theta <- theta
    z$lambda <- lambda
    
    errorhat <- rep(0,n) # E(error)=0 
    etahat <- rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i] <- beta0 + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]
     
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    del <- c(rep(1,n-m))
    mT <- diag(mu.eta(etahat[(m+1):n]))
    
    v <- ((z$lambda*exp(muhat^z$lambda)*muhat^(z$lambda-1)))/(1-exp(muhat^z$lambda))+
      ((z$lambda*log(0.5)*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*(1-exp(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2))
    
    c <- ((muhat^(z$lambda)*exp(muhat^z$lambda)*log(muhat))/(1-exp(muhat^z$lambda))) + (1/z$lambda)+
      ((log(0.5)*exp(muhat^z$lambda)*(muhat^z$lambda)*log(muhat)*(1-exp(y1^muhat)))/(1-exp(muhat^z$lambda))^2)+
      log(y1) - ((log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*log(y1))/(1-exp(muhat^z$lambda))) + ((y1^z$lambda)*log(y1))
    
    w <- (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*(z$lambda + z$lambda*(muhat^z$lambda) - (z$lambda-1)*exp(muhat^z$lambda) - 1))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda-1)*z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda^2)*exp(muhat^z$lambda)*(muhat^(2*z$lambda - 2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      (2*(z$lambda^2)*exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^3)
    
    d <- (exp(muhat^z$lambda)*(muhat^(z$lambda-1)) + z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat) + z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/(1 - exp(muhat^z$lambda))+
      (z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/((1 - exp(muhat^z$lambda))^2)-
      (exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*(muhat^(z$lambda-1))*log(0.5)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (2*z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^3)
    
    e <- -(log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*((y1^z$lambda)+1)*(log(y1)^2))/(1-exp(muhat^z$lambda)) + ((y1^z$lambda)*(log(y1)^2)) - (1/(z$lambda^2)) +
      (exp(muhat^z$lambda)*(muhat^z$lambda)*((muhat^z$lambda) + 1)*(log(muhat)^2))/(1-exp(muhat^z$lambda))+
      (exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^2)-
      (2*log(0.5)*(muhat^z$lambda)*log(muhat)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2)+
      (log(0.5)*(1-exp(y1^z$lambda))*exp(muhat^z$lambda)*(muhat^z$lambda)*(log(muhat)^2)*(1 + (muhat^z$lambda) ) )/((1-exp(muhat^z$lambda))^2)+
      (2*log(0.5)*(1-exp(y1^z$lambda))*exp(2*(muhat^lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^3)
    
    Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      Q[i,] <- errorhat[i+m-ma]
    }
    
    ###FB  recorrences
    deta.dbeta0 <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta <- matrix(0,ncol=q1,nrow=n)
    deta.dbeta0.dtheta <- matrix(0,ncol=q1,nrow=n)
    delta.theta.theta <- matrix(0,ncol=q1,nrow=q1)
    delta.theta.phi <- matrix(0,ncol=q1,nrow=p1)
    
    for(i in (m+1):n)
    {
      deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
      deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }
    
    for(j in 1:q1)
    {
      for(i in (m+1):n)
      {
        deta.dbeta0.dtheta[i,j] <- -deta.dbeta0[i-j] - theta%*%deta.dbeta0.dtheta[i-ma,j] 
      }
    }
    
    for(t in 1:q1){
      deta.dtheta.dtheta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dtheta[i,j] <- -deta.dtheta[i-t,j] - deta.dtheta[i-j, t] - theta%*%deta.dtheta.dtheta[i-ma,j]  
        }
      }   
      delta.theta.theta[t,] <- v%*%mT%*%deta.dtheta.dtheta[(m+1):n,]
    }
    
    for(t in 1:p1){
      deta.dphi.dtheta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dphi.dtheta[i,j] <- -deta.dphi[i-j, t] - theta%*%deta.dphi.dtheta[i-ma,j]  
        }
      }   
      delta.theta.phi[t,] <- v%*%mT%*%deta.dphi.dtheta[(m+1):n,]
    }
    
    b0 <- deta.dbeta0[(m+1):n]
    rQ <- deta.dtheta[(m+1):n,]
    rP <- deta.dphi[(m+1):n,]
    Aq <- deta.dbeta0.dtheta[(m+1):n,]
    
    W <- diag(w)
    V <- diag(v)
    D <- diag(d)
    E <- diag(e)
    
    ## Observed Information Matrix
    Kaa <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% b0 ##beta0_beta0
    Kap <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rP ##beta0_phi
    Kpa <- t(Kap) ##beta0_phi
    Kaq <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rQ -  v %*%mT%*%Aq ##beta0_theta
    Kqa <- t(Kaq) ##beta0_theta
    Kaprec <- - t(b0) %*% mT %*% D %*% del ##beta0_lambda
    Kpreca <- t(Kaprec) ##beta0_lambda
    Kpp <- -t(rP) %*% (W%*%(mT^2) + V%*%mT) %*% (rP) ##phi_phi
    Kpq <- -t(rP) %*% (W%*%(mT^2) + V%*%mT) %*% (rQ) - delta.theta.phi ##phi_theta
    Kqp <- t(Kpq) ##phi_theta
    Kpprec <- -t(rP) %*% mT %*% D %*% del ##phi_lambda
    Kprecp <- t(Kpprec) ##phi_lambda
    Kqq <- -t(rQ) %*% (W%*%(mT^2) + V%*%mT) %*% rQ - delta.theta.theta ##theta_theta
    Kqprec <- -t(rQ) %*% mT %*% D %*% del ##theta_lambda
    Kprecq <- t(Kqprec) ##theta_lambda
    Kprecprec <- -sum(e) ##lambda_lambda
    
    K <- rbind(
      cbind(Kaa,Kap,Kaq,Kaprec),
      cbind(Kpa,Kpp,Kpq,Kpprec),
      cbind(Kqa,Kqp,Kqq,Kqprec),
      cbind(Kpreca,Kprecp,Kprecq,Kprecprec)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    if(h1 != 0){
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- beta0 + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
        y_prev[n+i] <- linkinv(ynew_prev[n+i])
        errorhat[n+i] <- 0 
      }
    }
    
  }
  
  #   #####  AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T))
  {
    q1<-0
    reg <- c(mqo, lambda) # initializing the parameter values
    
    loglik <- function(z)
    {
      beta0 <- z[1]
      phi = z[2:(p1+1)]
      lambda <- z[p1+2] # precision parameter
      
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- beta0 + (phi%*%ynew[i-ar])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ETA <- (log(1- tau))/(1-(exp(mu^lambda))) 
      ll <- suppressWarnings( log( dchen(y1, b=lambda, lambda=ETA) ) )
      sum(ll)
    }
    
    escore <- function(z)
    {
      beta0 <- z[1]
      phi = z[2:(p1+1)]
      lambda <- z[p1+2] # precision parameter
      
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + (phi%*%ynew[i-ar])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      b0 <- c(rep(1,n-m))
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      v <- ((lambda*exp(mu^lambda)*mu^(lambda-1)))/(1-exp(mu^lambda))+
        ((lambda*log(0.5)*exp(mu^lambda)*(mu^(lambda-1))*(1-exp(y1^lambda)))/((1-exp(mu^lambda))^2))
      
      c <- ((mu^(lambda)*exp(mu^lambda)*log(mu))/(1-exp(mu^lambda))) + (1/lambda)+
        ((log(0.5)*exp(mu^lambda)*(mu^lambda)*log(mu)*(1-exp(y1^lambda)))/(1-exp(mu^lambda))^2)+
        log(y1) - ((log(0.5)*exp(y1^lambda)*(y1^lambda)*log(y1))/(1-exp(mu^lambda))) + ((y1^lambda)*log(y1))
      
      Ubeta0 <- t(b0) %*% mT %*% v
      Uphi <- t(P) %*% mT %*% v
      Ulambda <- sum(c)
      
      rval <- c(Ubeta0,Uphi,Ulambda)
      return(rval)
    }
    
    names_phi<-c(paste("phi",ar,sep=""))
    names_par <- c("beta0",names_phi,"lambda")
    
    opt <- optim(reg, loglik, gr=escore, hessian=T, method = "BFGS", 
                 control = list(fnscale = -1))
  
    z <- c()
    z$conv <- opt$conv
    coef <-c(opt$par)[1:(p1+2)]
    names(coef) <- names_par
    z$coeff <- coef
    
    beta0 <- coef[1]
    phi <- coef[2:(p1+1)]
    lambda <- coef[p1+2] # precision parameter
    
    z$beta0 <- beta0
    z$phi <- phi
    z$lambda <- lambda
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i] <- beta0 + (phi%*%ynew[i-ar])
      errorhat[i] <- ynew[i] - etahat[i] # predictor scale
    }
    
    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    b0 <- c(rep(1,n-m))
    del <- c(rep(1,n-m))
    mT <- diag(mu.eta(etahat[(m+1):n]))
    
    v <- ((z$lambda*exp(muhat^z$lambda)*muhat^(z$lambda-1)))/(1-exp(muhat^z$lambda))+
      ((z$lambda*log(0.5)*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*(1-exp(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2))
    
    c <- ((muhat^(z$lambda)*exp(muhat^z$lambda)*log(muhat))/(1-exp(muhat^z$lambda))) + (1/z$lambda)+
      ((log(0.5)*exp(muhat^z$lambda)*(muhat^z$lambda)*log(muhat)*(1-exp(y1^muhat)))/(1-exp(muhat^z$lambda))^2)+
      log(y1) - ((log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*log(y1))/(1-exp(muhat^z$lambda))) + ((y1^z$lambda)*log(y1))
    
    w <- (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*(z$lambda + z$lambda*(muhat^z$lambda) - (z$lambda-1)*exp(muhat^z$lambda) - 1))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda-1)*z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda^2)*exp(muhat^z$lambda)*(muhat^(2*z$lambda - 2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      (2*(z$lambda^2)*exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^3)
    
    d <- (exp(muhat^z$lambda)*(muhat^(z$lambda-1)) + z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat) + z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/(1 - exp(muhat^z$lambda))+
      (z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/((1 - exp(muhat^z$lambda))^2)-
      (exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*(muhat^(z$lambda-1))*log(0.5)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (2*z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^3)
    
    e <- -(log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*((y1^z$lambda)+1)*(log(y1)^2))/(1-exp(muhat^z$lambda)) + ((y1^z$lambda)*(log(y1)^2)) - (1/(z$lambda^2)) +
      (exp(muhat^z$lambda)*(muhat^z$lambda)*((muhat^z$lambda) + 1)*(log(muhat)^2))/(1-exp(muhat^z$lambda))+
      (exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^2)-
      (2*log(0.5)*(muhat^z$lambda)*log(muhat)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2)+
      (log(0.5)*(1-exp(y1^z$lambda))*exp(muhat^z$lambda)*(muhat^z$lambda)*(log(muhat)^2)*(1 + (muhat^z$lambda) ) )/((1-exp(muhat^z$lambda))^2)+
      (2*log(0.5)*(1-exp(y1^z$lambda))*exp(2*(muhat^lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^3)
    
    W <- diag(w)
    V <- diag(v)
    D <- diag(d)
    E <- diag(e)
    
    ## Observed Information Matrix
    Kaa <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% b0 ##beta0_beta0
    Kap <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% P ##beta0_phi
    Kpa <- t(Kap) ##beta0_phi
    Kaprec <- - t(b0) %*% mT %*% D %*% del ##beta0_lambda
    Kpreca <- t(Kaprec) ##beta0_lambda
    Kpp <- - t(P) %*% (W%*%(mT^2) + V%*%mT) %*% P ##phi_phi
    Kpprec <- -t(P) %*% mT %*% D %*% del ##phi_lambda
    Kprecp <- t(Kpprec) ##phi_lambda
    Kprecprec <- -sum(e) ##lambda_lambda
    
    K <- rbind(
      cbind(Kaa,Kap,Kaprec),
      cbind(Kpa,Kpp,Kpprec),
      cbind(Kpreca,Kprecp,Kprecprec)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    if(h1 != 0){
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- beta0 + (phi%*%ynew_prev[n+i-ar])
        y_prev[n+i] <- linkinv(ynew_prev[n+i])
        print(i)
      }
    }
  }
  
  #   ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T))
  {
    p1 <- 0
    reg <- c(mqo,rep(0,q1), lambda) # initializing the parameter values
    
    loglik <- function(z)
    {
      beta0 <- z[1]
      theta <- z[2:(q1+1)]
      lambda <- z[q1+2] # precision parameter
      
      eta <- error <- rep(0,n) # E(error)=0
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + (theta%*%error[i-ma])
        error[i] <- ynew[i] - eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
    
      ETA <- (log(1- tau))/(1-(exp(mu^lambda))) 
      ll <- suppressWarnings( log( dchen(y1, b=lambda, lambda=ETA) ) )
      sum(ll)
    }
    
    escore <- function(z)
    {
      beta0 <- z[1]
      theta <- z[2:(q1+1)]
      lambda <- z[q1+2] # precision parameter
      
      eta <- error <- rep(0,n) # E(error)=0
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + (theta%*%error[i-ma])
        error[i] <- ynew[i] - eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
      
      Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        Q[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrences
      deta.dbeta0 <- rep(0,n)
      deta.dtheta <- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
        deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }
      
      b0 <- deta.dbeta0[(m+1):n]
      qQ <- deta.dtheta[(m+1):n,]
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      v <- ((lambda*exp(mu^lambda)*mu^(lambda-1)))/(1-exp(mu^lambda))+
        ((lambda*log(0.5)*exp(mu^lambda)*(mu^(lambda-1))*(1-exp(y1^lambda)))/((1-exp(mu^lambda))^2))
      
      c <- ((mu^(lambda)*exp(mu^lambda)*log(mu))/(1-exp(mu^lambda))) + (1/lambda)+
        ((log(0.5)*exp(mu^lambda)*(mu^lambda)*log(mu)*(1-exp(y1^lambda)))/(1-exp(mu^lambda))^2)+
        log(y1) - ((log(0.5)*exp(y1^lambda)*(y1^lambda)*log(y1))/(1-exp(mu^lambda))) + ((y1^lambda)*log(y1))
      
      Ubeta0 <- t(b0) %*% mT %*% v
      Utheta <- t(qQ) %*% mT %*% v
      Ulambda <- sum(c)            
      
      rval <- c(Ubeta0,Utheta,Ulambda)
      return(rval)
    }
    names_theta<-c(paste("theta",ma,sep=""))
    names_par <- c("beta0",names_theta,"lambda")
    
    opt <- optim(reg, loglik, gr=escore, hessian=T, method = "BFGS", 
                 control = list(fnscale = -1))
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+2)]
    names(coef) <- names_par
    z$coeff <- coef
    
    beta0 <-coef[1]
    theta <- coef[2:(q1+1)]
    lambda <- coef[q1+2] # precision parameter
    
    z$beta0 <- beta0
    z$theta <- theta
    z$lambda <- lambda
    
    errorhat <- rep(0,n) # E(error)=0
    etahat <- rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i] <- beta0 + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i] - etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    del <- c(rep(1,n-m))
    mT <- diag(mu.eta(etahat[(m+1):n]))
    
    v <- ((z$lambda*exp(muhat^z$lambda)*muhat^(z$lambda-1)))/(1-exp(muhat^z$lambda))+
      ((z$lambda*log(0.5)*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*(1-exp(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2))
    
    c <- ((muhat^(z$lambda)*exp(muhat^z$lambda)*log(muhat))/(1-exp(muhat^z$lambda))) + (1/z$lambda)+
      ((log(0.5)*exp(muhat^z$lambda)*(muhat^z$lambda)*log(muhat)*(1-exp(y1^muhat)))/(1-exp(muhat^z$lambda))^2)+
      log(y1) - ((log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*log(y1))/(1-exp(muhat^z$lambda))) + ((y1^z$lambda)*log(y1))
    
    w <- (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*(z$lambda + z$lambda*(muhat^z$lambda) - (z$lambda-1)*exp(muhat^z$lambda) - 1))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda-1)*z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda^2)*exp(muhat^z$lambda)*(muhat^(2*z$lambda - 2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      (2*(z$lambda^2)*exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^3)
    
    d <- (exp(muhat^z$lambda)*(muhat^(z$lambda-1)) + z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat) + z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/(1 - exp(muhat^z$lambda))+
      (z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/((1 - exp(muhat^z$lambda))^2)-
      (exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*(muhat^(z$lambda-1))*log(0.5)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (2*z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^3)
    
    e <- -(log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*((y1^z$lambda)+1)*(log(y1)^2))/(1-exp(muhat^z$lambda)) + ((y1^z$lambda)*(log(y1)^2)) - (1/(z$lambda^2)) +
      (exp(muhat^z$lambda)*(muhat^z$lambda)*((muhat^z$lambda) + 1)*(log(muhat)^2))/(1-exp(muhat^z$lambda))+
      (exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^2)-
      (2*log(0.5)*(muhat^z$lambda)*log(muhat)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2)+
      (log(0.5)*(1-exp(y1^z$lambda))*exp(muhat^z$lambda)*(muhat^z$lambda)*(log(muhat)^2)*(1 + (muhat^z$lambda) ) )/((1-exp(muhat^z$lambda))^2)+
      (2*log(0.5)*(1-exp(y1^z$lambda))*exp(2*(muhat^lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^3)
    
    
    Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      Q[i,] <- errorhat[i+m-ma]
    }
    
    ###FB  recorrences
    deta.dbeta0 <- rep(0,n)
    deta.dtheta<- matrix(0,ncol=q1,nrow=n)
    deta.dbeta0.dtheta <- matrix(0,ncol=q1,nrow=n)
    delta.theta.theta <- matrix(0,ncol=q1,nrow=q1)
    
    for(i in (m+1):n)
    {
      deta.dbeta0[i]<- 1 - theta%*%deta.dbeta0[i-ma]
      deta.dtheta[i,]<- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }
    
    for(j in 1:q1)
    {
      for(i in (m+1):n)
      {
        deta.dbeta0.dtheta[i,j] <- - deta.dbeta0[i-j] - theta%*%deta.dbeta0.dtheta[i-ma,j] 
      }
    }
    
    for(t in 1:q1){
      deta.dtheta.dtheta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dtheta[i,j] <- -deta.dtheta[i-t,j] - deta.dtheta[i-j, t] - theta%*%deta.dtheta.dtheta[i-ma,j]  
        }
      }   
      delta.theta.theta[t,] <- v%*%mT%*%deta.dtheta.dtheta[(m+1):n,]
    }
    
    b0 <- deta.dbeta0[(m+1):n]
    qQ <- deta.dtheta[(m+1):n,]
    Aq <- deta.dbeta0.dtheta[(m+1):n,]
    
    W <- diag(w)
    V <- diag(v)
    D <- diag(d)
    E <- diag(e)
    
    ## Observed Information Matrix
    Kaa <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% b0 ##beta0_beta0
    Kaq <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% qQ -  v %*%mT%*%Aq##beta0_theta
    Kqa <- t(Kaq) ##beta0_phi
    Kaprec <- - t(b0) %*% mT %*% D %*% del ##beta0_lambda
    Kpreca <- t(Kaprec) ##beta0_lambda
    Kqq <- -t(qQ) %*% (W%*%(mT^2) + V%*%mT) %*% qQ - delta.theta.theta ##theta_theta
    Kqprec <- -t(qQ) %*% mT %*% D %*% del ##theta_lambda
    Kprecq <- t(Kqprec) ##phi_lambda
    Kprecprec <- -sum(e) ##lambda_lambda
    
    K <- rbind(
      cbind(Kaa,Kaq,Kaprec),
      cbind(Kqa,Kqq,Kqprec),
      cbind(Kpreca,Kprecq,Kprecprec)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    
    if(h1 != 0){
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- beta0 + (theta%*%errorhat[n+i-ma])
        y_prev[n+i] <- linkinv(ynew_prev[n+i])
        errorhat[n+i] <- 0 # original scale
      }
    }
    
    
  }
  
  ######### CHARMAX model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    beta1 <- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], rep(0,q1), lambda, beta1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      beta0 <- z[1]
      phi <- z[2:(p1+1)] 
      theta <- z[(p1+2):(p1+q1+1)]
      lambda <- z[p1+q1+2] # precision parameter
      beta <- z[(p1+q1+3):length(z)]
      
      error <- rep(0,n) # E(error)=0 
      eta <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] # predictor scale
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
      
      ETA <- (log(1- tau))/(1-(exp(mu^lambda))) 
      ll <- suppressWarnings(log(dchen(y1, b=lambda, lambda=ETA)))
      sum(ll)
    } 
    
    escore <- function(z)
    {
      beta0 <- z[1]
      phi <- z[2:(p1+1)] 
      theta <- z[(p1+2):(p1+q1+1)]
      lambda <- z[p1+q1+2] # precision parameter
      beta <- z[(p1+q1+3):length(z)]
      
      error <- rep(0,n) # E(error)=0 
      eta <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] # predictor scale
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n] 
      
      
      Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        Q[i,] <- error[i+m-ma]
      }
      
      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m))
      {
        P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
      }
    
      ###FB  recorrences
      deta.dbeta0 <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta <- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
        deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,] <- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      }
      
      b0 <- deta.dbeta0[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rQ <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      v <- ((lambda*exp(mu^lambda)*mu^(lambda-1)))/(1-exp(mu^lambda))+
        ((lambda*log(0.5)*exp(mu^lambda)*(mu^(lambda-1))*(1-exp(y1^lambda)))/((1-exp(mu^lambda))^2))
      
      c <- ((mu^(lambda)*exp(mu^lambda)*log(mu))/(1-exp(mu^lambda))) + (1/lambda)+
        ((log(0.5)*exp(mu^lambda)*(mu^lambda)*log(mu)*(1-exp(y1^lambda)))/(1-exp(mu^lambda))^2)+
        log(y1) - ((log(0.5)*exp(y1^lambda)*(y1^lambda)*log(y1))/(1-exp(mu^lambda))) + ((y1^lambda)*log(y1))
      
      Ubeta0 <- t(b0) %*% mT %*% v
      Uphi <-   t(rP) %*% mT %*% v
      Utheta <- t(rQ) %*% mT %*% v
      Ulambda <- sum(c)
      Ubeta <- t(rM) %*% mT %*% v
      
      rval <- c(Ubeta0,Uphi,Utheta,Ulambda,Ubeta)
      return(rval)
      
    }
    size.beta <- seq(1:length(X[1,]))
    names_beta <- c(paste("beta",size.beta,sep=""))
    names_phi <- c(paste("phi",ar,sep=""))
    names_theta <- c(paste("theta",ma,sep=""))
    names_par <- c("beta0",names_phi,names_theta,"lambda",names_beta)
    
    opt <- optim(reg, loglik, gr= escore, hessian=T, method = "BFGS", 
                 control = list(fnscale = -1))
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+2+ncol(X))]
    names(coef) <- names_par
    z$coeff <- coef
    
    beta0 <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    lambda <- coef[p1+q1+2] # precision parameter
    beta <- coef[(p1+q1+3):length(coef)]
    
    z$beta0 <- beta0
    z$phi <- phi
    z$theta <- theta
    z$lambda <- lambda
    z$beta <- beta
    
    errorhat <- rep(0,n) # E(error)=0 
    etahat <- rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i] <- beta0 + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i] - etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    ################################################
    del <- c(rep(1,n-m))
    mT <- diag(mu.eta(etahat[(m+1):n]))
    
    v <- ((z$lambda*exp(muhat^z$lambda)*muhat^(z$lambda-1)))/(1-exp(muhat^z$lambda))+
      ((z$lambda*log(0.5)*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*(1-exp(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2))
    
    c <- ((muhat^(z$lambda)*exp(muhat^z$lambda)*log(muhat))/(1-exp(muhat^z$lambda))) + (1/z$lambda)+
      ((log(0.5)*exp(muhat^z$lambda)*(muhat^z$lambda)*log(muhat)*(1-exp(y1^muhat)))/(1-exp(muhat^z$lambda))^2)+
      log(y1) - ((log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*log(y1))/(1-exp(muhat^z$lambda))) + ((y1^z$lambda)*log(y1))
    
    w <- (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*(z$lambda + z$lambda*(muhat^z$lambda) - (z$lambda-1)*exp(muhat^z$lambda) - 1))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda-1)*z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda^2)*exp(muhat^z$lambda)*(muhat^(2*z$lambda - 2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      (2*(z$lambda^2)*exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^3)
    
    d <- (exp(muhat^z$lambda)*(muhat^(z$lambda-1)) + z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat) + z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/(1 - exp(muhat^z$lambda))+
      (z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/((1 - exp(muhat^z$lambda))^2)-
      (exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*(muhat^(z$lambda-1))*log(0.5)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (2*z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^3)
    
    e <- -(log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*((y1^z$lambda)+1)*(log(y1)^2))/(1-exp(muhat^z$lambda)) + ((y1^z$lambda)*(log(y1)^2)) - (1/(z$lambda^2)) +
      (exp(muhat^z$lambda)*(muhat^z$lambda)*((muhat^z$lambda) + 1)*(log(muhat)^2))/(1-exp(muhat^z$lambda))+
      (exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^2)-
      (2*log(0.5)*(muhat^z$lambda)*log(muhat)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2)+
      (log(0.5)*(1-exp(y1^z$lambda))*exp(muhat^z$lambda)*(muhat^z$lambda)*(log(muhat)^2)*(1 + (muhat^z$lambda) ) )/((1-exp(muhat^z$lambda))^2)+
      (2*log(0.5)*(1-exp(y1^z$lambda))*exp(2*(muhat^lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^3)
    
    
    Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      Q[i,] <- errorhat[i+m-ma]
    }

    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }
    
    k1<- length(beta)
    M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
    for(i in 1:(n-m))
    {
      for(j in 1:length(beta))
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])
    }
    
    ###FB  recorrences
    deta.dbeta0 <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0,ncol=q1,nrow=n)
    deta.dbeta <- matrix(0, ncol=k1,nrow=n)
    deta.dbeta0.dtheta <- matrix(0,ncol=q1,nrow=n)
    delta.theta.theta <- matrix(0,ncol=q1,nrow=q1)
    delta.theta.phi <- matrix(0,ncol=q1,nrow=p1)
    delta.theta.beta <- matrix(0,ncol=q1,nrow=k1)
    delta.phi.beta <- matrix(0,ncol=p1,nrow=k1)
    
    for(i in (m+1):n)
    {
      deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
      deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,] <- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }
  
    for(j in 1:q1)
    {
      for(i in (m+1):n)
      {
        deta.dbeta0.dtheta[i,j] <- -deta.dbeta0[i-j] - theta%*%deta.dbeta0.dtheta[i-ma,j] 
      }
    }
    
    for(t in 1:q1){
      deta.dtheta.dtheta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dtheta[i,j] <- -deta.dtheta[i-t,j] - deta.dtheta[i-j, t] - theta%*%deta.dtheta.dtheta[i-ma,j]  
        }
      }   
      delta.theta.theta[t,] <- v%*%mT%*%deta.dtheta.dtheta[(m+1):n,]
    }
    
    for(t in 1:p1){
      deta.dphi.dtheta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dphi.dtheta[i,j] <- -deta.dphi[i-j, t] - theta%*%deta.dphi.dtheta[i-ma,j]  
        }
      }   
      delta.theta.phi[t,] <- v%*%mT%*%deta.dphi.dtheta[(m+1):n,]
    }
    
    for(t in 1:k1){
      deta.dtheta.dbeta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dbeta[i,j] <- -deta.dbeta[i-j, t] - theta%*%deta.dtheta.dbeta[i-ma,j]  
        }
      }   
      delta.theta.beta[t,] <- v%*%mT%*%deta.dtheta.dbeta[(m+1):n,]
    }
    
    for(t in 1:k1){
      deta.dbeta.dphi <- matrix(0,ncol=p1,nrow=n)
      for(j in 1:p1)
      {
        for(i in (m+1):n)
        {
          deta.dbeta.dphi[i,j] <- -(X[i-j,t]) - theta%*%deta.dbeta.dphi[i-ma,j]
        }
      }   
      delta.phi.beta[t,] <- v%*%mT%*%deta.dbeta.dphi[(m+1):n,]
    }
    
    b0 <- deta.dbeta0[(m+1):n]
    rQ <- deta.dtheta[(m+1):n,]
    rP <- deta.dphi[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    Aq <- deta.dbeta0.dtheta[(m+1):n,]
    
    W <- diag(w)
    V <- diag(v)
    D <- diag(d)
    E <- diag(e)
    
    ## Observed Information Matrix
    Kaa <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% b0 ##beta0_beta0
    Kap <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rP ##beta0_phi
    Kpa <- t(Kap) ##beta0_phi
    Kaq <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rQ -  v %*%mT%*%Aq ##beta0_theta
    Kqa <- t(Kaq) ##beta0_theta
    Kaprec <- -t(b0) %*% mT %*% D %*% del ##beta0_lambda
    Kpreca <- t(Kaprec) ##beta0_lambda
    Kab <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rM ##beta0_beta
    Kba <- t(Kab) ##beta0_beta
    Kpp <- -t(rP) %*% (W%*%(mT^2) + V%*%mT) %*% (rP) ##phi_phi
    Kpq <- -t(rP) %*% (W%*%(mT^2) + V%*%mT) %*% (rQ) - delta.theta.phi ##phi_theta
    Kqp <- t(Kpq) ##phi_theta
    Kpprec <- -t(rP) %*% mT %*% D %*% del ##phi_lambda
    Kprecp <- t(Kpprec) ##phi_lambda
    Kpb <- -t(rP) %*% (W%*%(mT^2) + V%*%mT) %*% (rM) - t(delta.phi.beta) ##phi_beta
    Kbp <- t(Kpb) ##phi_beta
    Kqq <- -t(rQ) %*% (W%*%(mT^2) + V%*%mT) %*% rQ - delta.theta.theta ##theta_theta
    Kqprec <- -t(rQ) %*% mT %*% D %*% del ##theta_lambda
    Kprecq <- t(Kqprec) ##theta_lambda
    Kqb <- -t(rQ) %*% (W%*%(mT^2) + V%*%mT) %*% rM - t(delta.theta.beta) ##theta_beta
    Kbq <- t(Kqb)
    Kprecprec <- -sum(e) ##lambda_lambda
    Kprecb <- -t(del) %*% mT %*% D %*% rM ##lambda_beta
    Kbprec <- t(Kprecb)
    Kbb <- -t(rM) %*% (W%*%(mT^2) + V%*%mT) %*% rM
    
    K <- rbind(
      cbind(Kaa,Kap,Kaq,Kaprec,Kab),
      cbind(Kpa,Kpp,Kpq,Kpprec,Kpb),
      cbind(Kqa,Kqp,Kqq,Kqprec,Kqb),
      cbind(Kpreca,Kprecp,Kprecq,Kprecprec,Kprecb),
      cbind(Kba,Kbp,Kbq,Kbprec,Kbb)
    )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    if(h1 != 0){
      X_prev<- rbind(X,X_hat)
      
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- beta0 + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
        y_prev[n+i] <- linkinv(ynew_prev[n+i])
        errorhat[n+i] <- 0 
      }
    }
    
     
  }
  
  ############# CHARX model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F))
  { 
    q1 <- 0
    beta1 <- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], lambda, beta1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      beta0 <- z[1]
      phi <- z[2:(p1+1)] 
      lambda <- z[p1+2] # precision parameter
      beta <- z[(p1+3):length(z)]
      
      error <- rep(0,n) # E(error)=0 
      eta <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
        error[i] <- ynew[i] - eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <-y[(m+1):n]
      
      ETA <- (log(1- tau))/(1-(exp(mu^lambda))) 
      ll <- suppressWarnings(log(dchen(y1, b=lambda, lambda=ETA)))
      sum(ll)
    }
    
    escore <- function(z)
    {
      beta0 <- z[1]
      phi <- z[2:(p1+1)] 
      lambda <- z[p1+2] # precision parameter
      beta <- z[(p1+3):length(z)]
      
      eta <- error <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
        error[i] <- ynew[i] - eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
      
      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m))
      {
        P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }
      
      M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          M[i,j] <- X[i+m,j] - sum(phi*X[i+m-ar,j])
      }
      
      b0 <- c(rep(1,n-m))
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      v <- ((lambda*exp(mu^lambda)*mu^(lambda-1)))/(1-exp(mu^lambda))+
        ((lambda*log(0.5)*exp(mu^lambda)*(mu^(lambda-1))*(1-exp(y1^lambda)))/((1-exp(mu^lambda))^2))
      
      c <- ((mu^(lambda)*exp(mu^lambda)*log(mu))/(1-exp(mu^lambda))) + (1/lambda)+
        ((log(0.5)*exp(mu^lambda)*(mu^lambda)*log(mu)*(1-exp(y1^lambda)))/(1-exp(mu^lambda))^2)+
        log(y1) - ((log(0.5)*exp(y1^lambda)*(y1^lambda)*log(y1))/(1-exp(mu^lambda))) + ((y1^lambda)*log(y1))
      
      Ubeta0 <- t(b0) %*% mT %*% v
      Ubeta <- t(M) %*% mT %*% v
      Uphi <- t(P) %*% mT %*% v
      Ulambda <- sum(c)
      
      rval <- c(Ubeta0,Uphi,Ulambda,Ubeta)
      return(rval)
    }
    size.beta <- seq(1:length(X[1,])) 
    names_phi <- c(paste("phi",ar,sep=""))
    names_beta <- c(paste("beta",size.beta,sep=""))
    names_par <- c("beta0",names_phi,"lambda",names_beta)
    
    opt <- optim(reg, loglik, gr=escore, hessian=T, method = "BFGS", 
                 control = list(fnscale = -1))
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+2+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    beta0 <-coef[1]
    phi <- coef[2:(p1+1)]
    lambda <- coef[p1+2] # precision parameter
    beta <- coef[(p1+3):length(coef)]
    
    z$beta0 <- beta0
    z$phi <- phi
    z$lambda <- lambda
    z$beta <- beta0
    
    errorhat <- rep(0,n) # E(error)=0 
    etahat <- rep(NA,n)
    
    for(i in (m+1):n){
      etahat[i] <- beta0 + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
      errorhat[i] <- ynew[i] - etahat[i] # predictor scale
    }
    
    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]
     
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }
    
    M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
    for(i in 1:(n-m))
    {
      for(j in 1:length(beta))
        M[i,j] <- X[i+m,j] - sum(phi*X[i+m-ar,j])
    }
     
    b0 <- c(rep(1,n-m))
    del <- c(rep(1,n-m))
    mT <- diag(mu.eta(etahat[(m+1):n]))
     
    v <- ((z$lambda*exp(muhat^z$lambda)*muhat^(z$lambda-1)))/(1-exp(muhat^z$lambda))+
       ((z$lambda*log(0.5)*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*(1-exp(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2))
     
    c <- ((muhat^(z$lambda)*exp(muhat^z$lambda)*log(muhat))/(1-exp(muhat^z$lambda))) + (1/z$lambda)+
       ((log(0.5)*exp(muhat^z$lambda)*(muhat^z$lambda)*log(muhat)*(1-exp(y1^muhat)))/(1-exp(muhat^z$lambda))^2)+
       log(y1) - ((log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*log(y1))/(1-exp(muhat^z$lambda))) + ((y1^z$lambda)*log(y1))
     
    w <- (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*(z$lambda + z$lambda*(muhat^z$lambda) - (z$lambda-1)*exp(muhat^z$lambda) - 1))/((1 - exp(muhat^z$lambda))^2)+
       ((z$lambda-1)*z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
       ((z$lambda^2)*exp(muhat^z$lambda)*(muhat^(2*z$lambda - 2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
       (2*(z$lambda^2)*exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^3)
     
    d <- (exp(muhat^z$lambda)*(muhat^(z$lambda-1)) + z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat) + z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/(1 - exp(muhat^z$lambda))+
       (z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/((1 - exp(muhat^z$lambda))^2)-
       (exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
       (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
       (z$lambda*(muhat^(z$lambda-1))*log(0.5)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)-
       (z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
       (2*z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^3)
     
    e <- -(log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*((y1^z$lambda)+1)*(log(y1)^2))/(1-exp(muhat^z$lambda)) + ((y1^z$lambda)*(log(y1)^2)) - (1/(z$lambda^2)) +
       (exp(muhat^z$lambda)*(muhat^z$lambda)*((muhat^z$lambda) + 1)*(log(muhat)^2))/(1-exp(muhat^z$lambda))+
       (exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^2)-
       (2*log(0.5)*(muhat^z$lambda)*log(muhat)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2)+
       (log(0.5)*(1-exp(y1^z$lambda))*exp(muhat^z$lambda)*(muhat^z$lambda)*(log(muhat)^2)*(1 + (muhat^z$lambda) ) )/((1-exp(muhat^z$lambda))^2)+
       (2*log(0.5)*(1-exp(y1^z$lambda))*exp(2*(muhat^lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^3)
     
    W <- diag(w)
    V <- diag(v)
    D <- diag(d)
    E <- diag(e)
     
    ## Observed Information Matrix
    Kaa <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% b0 ##beta0_beta0
    Kap <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% P ##beta0_phi
    Kpa <- t(Kap) ##beta0_phi
    Kaprec <- -t(b0) %*% mT %*% D %*% del ##beta0_lambda
    Kpreca <- t(Kaprec) ##beta0_lambda
    Kab <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% M ##beta0_beta
    Kba <- t(Kab) ##beta0_beta
    Kpp <- -t(P) %*% (W%*%(mT^2) + V%*%mT) %*% P ##phi_phi
    Kpprec <- -t(P) %*% mT %*% D %*% del ##phi_lambda
    Kprecp <- t(Kpprec) ##phi_lambda
    Kpb <- -t(P) %*% (W%*%(mT^2) + V%*%mT) %*% M ##phi_beta
    Kbp <- t(Kpb)
    Kprecprec <- -sum(e) ##lambda_lambda
    Kbprec <- -t(M) %*% mT %*% D %*% del ##beta0_lambda
    Kprecb <- t(Kbprec) ##beta0_lambda
    Kbb <- -t(M) %*% (W%*%(mT^2) + V%*%mT) %*% M
    
    K <- rbind(
      cbind(Kaa,Kap,Kaprec,Kab),
      cbind(Kpa,Kpp,Kpprec,Kpb),
      cbind(Kpreca,Kprecp,Kprecprec,Kprecb),
      cbind(Kba,Kbp,Kbprec,Kbb)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    if(h1 != 0){
      X_prev<- rbind(X,X_hat)
      
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- beta0 + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
        y_prev[n+i] <- linkinv(ynew_prev[n+i])
        errorhat[n+i] <- 0 
      }
    }
    
  }
  
  #   ######### CMAX model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F))
  {
    p1 <- 0
    beta1 <- mqo[(2):length(mqo)]
    reg <- c(mqo[1], rep(0,q1), lambda, beta1) # initializing the parameter values
    
    loglik <- function(z)
    {
      beta0 <- z[1]
      theta <- z[(2):(q1+1)]
      lambda <- z[q1+2] # precision parameter
      beta <- z[(q1+3):length(z)]
      
      error <- rep(0,n) # E(error)=0
      eta <- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
      
      ETA <- (log(1- tau))/(1-(exp(mu^lambda))) 
      ll <- suppressWarnings(log(dchen(y1, b=lambda, lambda=ETA)))
      sum(ll)
    }
    
    escore <- function(z)
    {
      beta0 <- z[1]
      theta <- z[(2):(q1+1)]
      lambda <- z[q1+2] # precision parameter
      beta <- z[(q1+3):length(z)]
      
      eta <- error <- rep(0,n) # E(error)=0
      
      for(i in (m+1):n)
      {
        eta[i] <- beta0 + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1 <- y[(m+1):n]
      
      Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        Q[i,] <- error[i+m-ma]
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+m,j]
      }
      
      ###FB  recorrencias
      deta.dbeta0 <- rep(0,n)
      deta.dtheta <- matrix(0, ncol=q1,nrow=n)
      deta.dbeta <- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
        deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,] <- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
        
      }
      
      b0 <- deta.dbeta0[(m+1):n]
      rQ <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      v <- ((lambda*exp(mu^lambda)*mu^(lambda-1)))/(1-exp(mu^lambda))+
        ((lambda*log(0.5)*exp(mu^lambda)*(mu^(lambda-1))*(1-exp(y1^lambda)))/((1-exp(mu^lambda))^2))
      
      c <- ((mu^(lambda)*exp(mu^lambda)*log(mu))/(1-exp(mu^lambda))) + (1/lambda)+
        ((log(0.5)*exp(mu^lambda)*(mu^lambda)*log(mu)*(1-exp(y1^lambda)))/(1-exp(mu^lambda))^2)+
        log(y1) - ((log(0.5)*exp(y1^lambda)*(y1^lambda)*log(y1))/(1-exp(mu^lambda))) + ((y1^lambda)*log(y1))
      
      Ubeta0 <- t(b0) %*% mT %*% v
      Utheta <- t(rQ) %*% mT %*% v
      Ulambda <- sum(c)    
      Ubeta <- t(rM) %*% mT %*% v
      
      rval <- c(Ubeta0,Utheta,Ulambda,Ubeta)
      return(rval)
    }
    size.beta <- seq(1:length(X[1,])) 
    names_beta <- c(paste("beta",size.beta,sep=""))
    names_theta<-c(paste("theta",ma,sep=""))
    names_par <- c("beta0",names_theta,"lambda",names_beta)
    
    opt <- optim(reg, loglik, gr=escore, hessian=T, method = "BFGS", 
                 control = list(fnscale = -1))
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+2+ncol(X) )]
    names(coef) <- names_par
    z$coeff <- coef
    
    beta0 <-coef[1]
    theta <- coef[(2):(q1+1)]
    lambda <- coef[q1+2] # precision parameter
    beta <- coef[(q1+3):length(coef)]
    
    z$beta0 <- beta0
    z$theta <- theta
    z$lambda <- lambda
    
    errorhat <- rep(0,n) # E(error)=0
    etahat <- rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i] <- beta0 + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    del <- c(rep(1,n-m))
    mT <- diag(mu.eta(etahat[(m+1):n]))
    
    v <- ((z$lambda*exp(muhat^z$lambda)*muhat^(z$lambda-1)))/(1-exp(muhat^z$lambda))+
      ((z$lambda*log(0.5)*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*(1-exp(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2))
    
    c <- ((muhat^(z$lambda)*exp(muhat^z$lambda)*log(muhat))/(1-exp(muhat^z$lambda))) + (1/z$lambda)+
      ((log(0.5)*exp(muhat^z$lambda)*(muhat^z$lambda)*log(muhat)*(1-exp(y1^muhat)))/(1-exp(muhat^z$lambda))^2)+
      log(y1) - ((log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*log(y1))/(1-exp(muhat^z$lambda))) + ((y1^z$lambda)*log(y1))
    
    w <- (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*(z$lambda + z$lambda*(muhat^z$lambda) - (z$lambda-1)*exp(muhat^z$lambda) - 1))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda-1)*z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      ((z$lambda^2)*exp(muhat^z$lambda)*(muhat^(2*z$lambda - 2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)+
      (2*(z$lambda^2)*exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda-2))*log(0.5)*(1-exp(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^3)
    
    d <- (exp(muhat^z$lambda)*(muhat^(z$lambda-1)) + z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat) + z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/(1 - exp(muhat^z$lambda))+
      (z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat))/((1 - exp(muhat^z$lambda))^2)-
      (exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*(muhat^(z$lambda-1))*log(0.5)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1 - exp(muhat^z$lambda))^2)-
      (z$lambda*exp(muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^2)-
      (2*z$lambda*exp(2*muhat^z$lambda)*(muhat^(2*z$lambda-1))*log(muhat)*log(0.5)*(exp(y1^z$lambda)-1))/((1 - exp(muhat^z$lambda))^3)
    
    e <- -(log(0.5)*exp(y1^z$lambda)*(y1^z$lambda)*((y1^z$lambda)+1)*(log(y1)^2))/(1-exp(muhat^z$lambda)) + ((y1^z$lambda)*(log(y1)^2)) - (1/(z$lambda^2)) +
      (exp(muhat^z$lambda)*(muhat^z$lambda)*((muhat^z$lambda) + 1)*(log(muhat)^2))/(1-exp(muhat^z$lambda))+
      (exp(2*(muhat^z$lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^2)-
      (2*log(0.5)*(muhat^z$lambda)*log(muhat)*(y1^z$lambda)*log(y1)*exp((muhat^z$lambda)+(y1^z$lambda)))/((1-exp(muhat^z$lambda))^2)+
      (log(0.5)*(1-exp(y1^z$lambda))*exp(muhat^z$lambda)*(muhat^z$lambda)*(log(muhat)^2)*(1 + (muhat^z$lambda) ) )/((1-exp(muhat^z$lambda))^2)+
      (2*log(0.5)*(1-exp(y1^z$lambda))*exp(2*(muhat^lambda))*(muhat^(2*z$lambda))*(log(muhat)^2))/((1-exp(muhat^z$lambda))^3)
    
    Q <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      Q[i,] <- errorhat[i+m-ma]
    }
    
    k1<-length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]
    }
    
    ###FB  recorrencias
    deta.dbeta0 <- rep(0,n)
    deta.dtheta <- matrix(0, ncol=q1,nrow=n)
    deta.dbeta <- matrix(0, ncol=k1,nrow=n)
    deta.dtheta.dbeta0 <- matrix(0,ncol=q1,nrow=n)
    delta.theta.theta <- matrix(0,ncol=q1,nrow=q1)
    delta.theta.beta <- matrix(0,ncol=q1,nrow=k1)
    
    for(i in (m+1):n)
    {
      deta.dbeta0[i] <- 1 - theta%*%deta.dbeta0[i-ma]
      deta.dtheta[i,] <- Q[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,] <- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      
    }
    
    for(j in 1:q1)
    {
      for(i in (m+1):n)
      {
        deta.dtheta.dbeta0[i,j] <- -deta.dbeta0[i-j] - theta%*%deta.dtheta.dbeta0[i-ma,j] 
      }
    }
    
    for(t in 1:q1){
      deta.dtheta.dtheta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dtheta[i,j] <- -deta.dtheta[i-t,j] - deta.dtheta[i-j, t] - theta%*%deta.dtheta.dtheta[i-ma,j]  
        }
      }   
      delta.theta.theta[t,] <- v%*%mT%*%deta.dtheta.dtheta[(m+1):n,]
    }
    
    
    for(t in 1:k1){
      deta.dtheta.dbeta <- matrix(0,ncol=q1,nrow=n)
      for(j in 1:q1)
      {
        for(i in (m+1):n)
        {
          deta.dtheta.dbeta[i,j] <- -deta.dbeta[i-j, t] - theta%*%deta.dtheta.dbeta[i-ma,j]  
        }
      }   
      delta.theta.beta[t,] <- v%*%mT%*%deta.dtheta.dbeta[(m+1):n,]
    }
    
    b0 <- deta.dbeta0[(m+1):n]
    Aq <- deta.dtheta.dbeta0[(m+1):n,]
    rQ <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
    W <- diag(w)
    V <- diag(v)
    D <- diag(d)
    E <- diag(e)
    
    ## Observed Information Matrix
    Kaa <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% b0 ##beta0_beta0
    Kaq <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rQ -  v %*%mT%*%Aq##beta0_theta
    Kqa <- t(Kaq) ##beta0_phi
    Kaprec <- -t(b0) %*% mT %*% D %*% del ##beta0_lambda
    Kpreca <- t(Kaprec) ##beta0_lambda
    Kab <- -t(b0) %*% (W%*%(mT^2) + V%*%mT) %*% rM ##beta0_beta
    Kba <- t(Kab) ##beta0_beta
    Kqq <- -t(rQ) %*% (W%*%(mT^2) + V%*%mT) %*% rQ - t(delta.theta.theta) ##theta_theta
    Kqprec <- -t(rQ) %*% mT %*% D %*% del ##theta_lambda
    Kprecq <- t(Kqprec) ##theta_lambda
    Kqb <- -t(rQ) %*% (W%*%(mT^2) + V%*%mT) %*% rM - t(delta.theta.beta) ##theta_theta
    Kbq <- t(Kqb)
    Kprecprec <- -sum(e) ##lambda_lambda
    Kprecb <- -t(del) %*% mT %*% D %*% rM ##lambda_beta
    Kbprec <- t(Kprecb)
    Kbb <- -t(rM) %*% (W%*%(mT^2) + V%*%mT) %*% rM
    
    K <- rbind(
      cbind(Kaa,Kaq,Kaprec,Kab),
      cbind(Kqa,Kqq,Kqprec,Kqb),
      cbind(Kpreca,Kprecq,Kprecprec,Kprecb),
      cbind(Kba,Kbq,Kbprec,Kbb)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    
    if(h1 != 0){
      X_prev<- rbind(X,X_hat)
      
      for(i in 1:h1)
      {
        ynew_prev[n+i] <- beta0 + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
        y_prev[n+i] <- linkinv(ynew_prev[n+i])
        errorhat[n+i] <- 0 
      }
    }
     
  }

  z$serie <- y
  z$barma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
   
  ETA <- (log(1- tau))/(1-(exp(muhat^lambda)))
  z$residual = qnorm(pchen(y[(m+1):n], b=lambda, lambda=ETA, log.p=FALSE, lower.tail=TRUE))
  resid <- z$residual
  
###########################################################################################  

   vcov <- chol2inv(chol(K))
   z$vcov <- vcov
   
   stderror <- sqrt(diag(vcov))
   z$stderror <- stderror
   
   z$zstat <- abs(z$coef/stderror)
   z$pvalues <- 2*(1 - pnorm(z$zstat) )
   
   z$loglik <- opt$value
   z$counts <- as.numeric(opt$counts[1])
   
  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+2+length(beta))
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }
   
   model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),
                                round(z$zstat,4),round(z$pvalues,4))
   colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
   z$model <- model_presentation
   
  if( desc==1){
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)),quote=F)
    print("Residuals:",quote=F)
    print(summary(as.vector(resid)))
    
  }
    
    ###################################################
    ######### GRAPHICS ################################
    
  if(graf>0)
    {
    # print(model_presentation)
    # print(" ",quote=F)
    # print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    # print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    # print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
    # print("Residuals:",quote=F)
    # print(summary(resid)) 
      
    t<-seq(-5,n+6,by=1)
      
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(resid,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)
      
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="Fitted values",ylab="Observed data",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
      
    plot(as.vector(z$fitted[(m+1):n]),as.vector(resid), main=" ", pch = "+",
         xlab="Fitted values",ylab="Residuals")
      #lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    densidade<-density(resid)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
      
    acf(resid,ylab="ACF",xlab="Lag") 
    pacf(resid,ylab="PACF",xlab="Lag") 
      
    max_r<- max(resid,na.rm=T)
    min_r<- min(resid,na.rm=T)
    qqnorm(resid, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="Normal quantiles",
           ylab="Empirical quantiles")
    lines(c(-10,10),c(-10,10),lty=2)
      
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="Serie",xlab="Time")
    lines(z$fitted,col="red")
      
    fim<-end(y)[1]+end(y)[2]/12
      
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    par(mfrow=c(1,1))
    plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
    abline(v=fim,lty=2)
    lines(y)
      
    w1<-5
    h1<-4
    
    if(graf>1)
      {
      postscript(file = "resid_v_ind.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(resid,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
        }
        dev.off()
        
      pdf(file = "resid_v_fitted.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted[(m+1):n]),as.vector(resid), main=" ", pch = "+",
             xlab="Fitted values",ylab="Residuals",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
        
      pdf(file = "obs_v_fit.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
             xlab="Fitted values",ylab="Observed data",
             xlim=c(0.95*min_y,max_y*1.05),
             ylim=c(0.95*min_y,max_y*1.05))
        lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
        }
        dev.off()
        
      pdf(file = "resid_density.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(1.5, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
          
        plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
        lines(densidade$x,dnorm(densidade$x),lty=2)
        legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n")
      }
      dev.off()
        
      pdf(file = "resid_FAC.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        acf(resid,ylab="ACF",xlab="Lag") 
      }
      dev.off()
        
      pdf(file = "resid_FACP.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        pacf(resid,ylab="PACF",xlab="Lag")
      }
      dev.off()
        
      pdf(file = "qq_plot.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(resid, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="Normal quantiles",ylab="Empirical quantiles")
        lines(c(-10,10),c(-10,10),lty=2)
      }
      dev.off()
        
      pdf(file = "adjusted.pdf",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="Serie",xlab="Time")
        lines(z$fitted,col="red")
      }
      dev.off()
        
      pdf(file = "forecast.pdf",horizontal=F,paper="special",width = 6, height = 4.7,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",lty=2,col="red", ylim=c(min(y),max(y)),ylab="RH",xlab="Times")
        abline(v=fim,lty=2)
        lines(y)
        legend("bottomleft",c("Observed data","Fitted and forecast values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,"red"))
      }
      dev.off()
        
      }    
    }  
    
    return(z)
}

  
  
  

