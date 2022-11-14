require(flexsurv)
require(numDeriv)

OLLMax.mixtura <- function (mu.link = "log", sigma.link="log", nu.link = "logit")
{
  mstats <- checklink(   "mu.link", "OLLMax.mixtura", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "OLLMax.mixtura", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink(   "nu.link", "OLLMax.mixtura", substitute(nu.link),    
                         c("1/nu^2", "log", "identity", "own"))
   structure(
    list(family = c("OLLMax.mixtura", "OLLMax.mixtura"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         
         dldm = function(y,mu,sigma,nu){ #----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLLMax.mixtura(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLLMax.mixtura(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         dldd = function(y,mu,sigma,nu){#----------------------------------------------------- ok  
           lpdf<-function(t,mu,x,nu){log(dauxiOLLMax.mixtura(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu){log(dauxiOLLMax.mixtura(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },   
         dldv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,x){log(dauxiOLLMax.mixtura(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,method='simple')
           dldv
         },
         d2ldv2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok 
           lpdf<-function(t,mu,sigma,x){log(dauxiOLLMax.mixtura(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,method='simple')
           d2ldv2<- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
           d2ldv2
         },
  
         d2ldmdd = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLLMax.mixtura(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,x,nu){log(dauxiOLLMax.mixtura(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },
         d2ldmdv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu){log(dauxiOLLMax.mixtura(t,x,sigma,nu))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,sigma,x){log(dauxiOLLMax.mixtura(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu)
           d2ldmdv = -(dldm * dldv)
           d2ldmdv				
         },
         
         d2ldddv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu){log(dauxiOLLMax.mixtura(t,mu,x,nu))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,method='simple')
           lpdf<-function(t,mu,sigma,x){log(dauxiOLLMax.mixtura(t,mu,sigma,x))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu)
           d2ldddv = -(dldd * dldv)
           d2ldddv	
         },
         #----------------------------------------------------- ok
         G.dev.incr  = function(y,mu,sigma,nu,...) 
         { 
           -2*dOLLMax.mixtura(y,mu,sigma,nu,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pOLLMax.mixtura", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)) ,
         #mu.initial = expression(mu <- (y+mean(y))/2), 
         mu.initial = expression(mu <- y+mean(y)),
         #mu.initial = expression(mu <- rep(2, length(y))),
         #mu.initial = expression(mu <- rep(sd(y), length(y))), 
         sigma.initial = expression(sigma <- rep(2.45, length(y))), #OK
         nu.initial = expression(nu <- rep(0.65,length(y))), #)k
         mu.valid = function(mu) all(mu > 0), 
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu >= 0.01) && all(nu <= .99),
         y.valid = function(y)  all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dOLLMax.mixtura  <- function(x, mu = 2, sigma = 1, nu = 0.5,log = FALSE){
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  
  
  g <-  dgengamma.orig(x,shape=2,scale=mu,k=3/2)
  G <- pgengamma.orig(x,shape=2,scale=mu,k=3/2)
  fy1 <- (1-nu)*(sigma*g*(G^(sigma-1))*(1-G)^(sigma-1))/(G^(sigma)+(1-G)^sigma)^2
  
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}    
#-----------------------------------------------------------------  
pOLLMax.mixtura <- function(q, mu = 2, sigma = 1, nu = 0.5, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  
  #if (any((tau >= 0.01) && all(tau < .99)))  stop(paste("tau must be positive", "\n", ""))
  
  G <- pgengamma.orig(q,shape=2,scale=mu,k=3/2)
  FF1 <- (G^(sigma))/((G^(sigma)+(1-G)^sigma))
  cdf1 <- (1-nu)*FF1
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}



#-----------------------------------------------------------------  

qOLLMax.mixtura <- function(q, mu = 2, sigma = 1, nu = 0.5, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  l1 <- 1/sigma
  u1 <- (p^(l1))/(((1-p)^l1)+p^l1)
  q <- qgengamma.orig(u1,scale =mu ,shape=2,k=3/2)
  q
}



#-----------------------------------------------------------------  


dauxiOLLMax.mixtura <- function(t,mu,sigma,nu){ 
  g <-  dgengamma.orig(t,shape=2,scale=mu,k=3/2)
  G <- pgengamma.orig(t,shape=2,scale=mu,k=3/2)
  fy1 <- (1-nu)*(sigma*g*(G^(sigma-1))*(1-G)^(sigma-1))/(G^(sigma)+(1-G)^sigma)^2
  fy1}

