require(flexsurv)
require(numDeriv)

GOLLMax.mixtura <- function (mu.link = "log", sigma.link="log", nu.link = "log", tau.link = "logit")
{
  mstats <- checklink(   "mu.link", "GOLLMax.mixtura", substitute(mu.link), 
                         c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "GOLLMax.mixtura", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink(   "nu.link", "GOLLMax.mixtura", substitute(nu.link),    
                         c("1/nu^2", "log", "identity", "own"))
  tstats <- checklink(  "tau.link", "GOLLMax.mixtura", substitute(tau.link),   
                        c("logit", "cauchit", "cloglog", "probit","loglog")) 
  structure(
    list(family = c("GOLLMax.mixtura", "GOLLMax.mixtura"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
         nopar = 4, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         tau.link = as.character(substitute(tau.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         tau.linkfun = tstats$linkfun,  
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         tau.linkinv = tstats$linkinv, 
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         tau.dr = tstats$mu.eta, 
         
         dldm = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiGOLLMax.mixtura(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiGOLLMax.mixtura(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           d2ldm2 <- -dldm * dldm
           d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
           d2ldm2 
         },     
         dldd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok  
           lpdf<-function(t,mu,x,nu,tau){log(dauxiGOLLMax.mixtura(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu,tau){log(dauxiGOLLMax.mixtura(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           d2ldd2 <- -dldd*dldd
           d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
           d2ldd2
         },   
         dldv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiGOLLMax.mixtura(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
           dldv
         },
         d2ldv2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok 
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiGOLLMax.mixtura(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau,method='simple')
           d2ldv2<- -dldv * dldv
           d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                   
           d2ldv2
         },
         dldt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiGOLLMax.mixtura(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           dldt
         } ,
         d2ldt2 = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiGOLLMax.mixtura(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau,method='simple')
           d2ldt2<- -dldt * dldt
           d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
           d2ldt2
         },
         d2ldmdd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiGOLLMax.mixtura(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,x,nu,tau){log(dauxiGOLLMax.mixtura(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
           d2ldmdd = -(dldm * dldd)
           d2ldmdd<-ifelse(d2ldmdd < -1e-15, d2ldmdd,-1e-15)
           #d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
           d2ldmdd                 
         },
         d2ldmdv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiGOLLMax.mixtura(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiGOLLMax.mixtura(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           d2ldmdv = -(dldm * dldv)
           d2ldmdv				
         },
         
         d2ldmdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,x,sigma,nu,tau){log(dauxiGOLLMax.mixtura(t,x,sigma,nu,tau))}
           dldm<-grad(func=lpdf,t=y,x=mu,sigma=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiGOLLMax.mixtura(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldmdt <- -(dldm*dldt)
           d2ldmdt
         },
         
         d2ldddv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu,tau){log(dauxiGOLLMax.mixtura(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau,method='simple')
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiGOLLMax.mixtura(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           d2ldddv = -(dldd * dldv)
           d2ldddv	
         },
         d2ldddt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,mu,x,nu,tau){log(dauxiGOLLMax.mixtura(t,mu,x,nu,tau))}
           dldd<-grad(func=lpdf,t=y,mu=mu,x=sigma,nu=nu,tau=tau)
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiGOLLMax.mixtura(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldddt <- -(dldd*dldt) 
           d2ldddt 
         },
         d2ldvdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           lpdf<-function(t,mu,sigma,x,tau){log(dauxiGOLLMax.mixtura(t,mu,sigma,x,tau))}
           dldv<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,x=nu,tau=tau)
           lpdf<-function(t,mu,sigma,nu,x){log(dauxiGOLLMax.mixtura(t,mu,sigma,nu,x))}
           dldt<-grad(func=lpdf,t=y,mu=mu,sigma=sigma,nu=nu,x=tau)
           d2ldvdt <- -(dldv*dldt) 
           d2ldvdt 
         },
         #----------------------------------------------------- ok
         G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
         { 
           -2*dGOLLMax.mixtura(y,mu,sigma,nu,tau,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pGOLLMax.mixtura", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
        #mu.initial = expression(mu <- rep(sd(y)*2.2+sqrt(pi), length(y))), 
        mu.initial = expression(mu <- y+mean(y)),
         sigma.initial = expression(sigma <- rep(1.2, length(y))), #OK
         nu.initial = expression(nu <- rep(2.45,length(y))), #)k
         tau.initial = expression(tau <-rep(0.65, length(y))), 
         mu.valid = function(mu) all(mu > 0), 
         sigma.valid = function(sigma) all(sigma > 0),
         nu.valid = function(nu) all(nu > 0) ,
         tau.valid = function(tau) all(tau >= 0.01) && all(tau <= .99),
         y.valid = function(y)  all(y > 0)
    ),
    class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dGOLLMax.mixtura  <- function(x, mu = 2, sigma = 1, nu = 0.5, tau = 0.5, log = FALSE){
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
 #if (any((tau >= 0.01) && all(tau < .99)))  stop(paste("tau must be between zero and one", "\n", ""))
  
  g <-  dgengamma.orig(x,shape=2,scale=mu,k=3/2)
  G <- pgengamma.orig(x,shape=2,scale=mu,k=3/2)
  f <- (1-tau)*(nu*sigma*g*(G^(nu*sigma-1))*(1-G^sigma)^(nu-1))/(G^(nu*sigma)+(1-G^sigma)^nu)^2
  
  if(log==FALSE) fy<-f else fy<-log(f)
  fy
}    
#-----------------------------------------------------------------  
pGOLLMax.mixtura <- function(q, mu = 2, sigma = 1, nu = 0.5, tau = 0.5, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", "")) 
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))  
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  #if (any((tau >= 0.01) && all(tau < .99)))  stop(paste("tau must be between zero and one", "\n", ""))
  
  G <- pgengamma.orig(q,shape=2,scale=mu,k=3/2)
  FF1 <- (G^(nu*sigma))/((G^(nu*sigma)+(1-G^sigma)^nu))
  cdf1 <- (1-tau)*FF1
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#-----------------------------------------------------------------  
#-----------------------------------------------------------------  

sGOLLMax.mixtura <- function(x, mu = 2, sigma = 1, nu = 0.5, tau = .5, lower.tail = TRUE, log.p = FALSE){  
 S=(1-pGOLLMax.mixtura(x,mu,sigma,nu,tau))
}

#-----------------------------------------------------------------  
hGOLLMax.mixtura <- function(x, mu = 2, sigma = 1, nu = 0.5, tau = .5, lower.tail = TRUE, log.p = FALSE){  
  hazard=dGOLLMax.mixtura(x,mu,sigma,nu,tau)/sGOLLMax.mixtura(x,mu,sigma,nu,tau)
}
#-----------------------------------------------------------------  
qGOLLMax.mixtura <-  function(p, mu = 2, sigma = 1, nu = 0.5, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  
  l1 <- 1/nu
  l2 <- 1/sigma
  u1 <- (p^(l1*l2))/(((1-p)^l1)+p^l1)^l2
  q <- qgengamma.orig(u1,scale =mu ,shape=2,k=3/2)
  q
}
#-----------------------------------------------------------------  
rGOLLMax.mixtura <- function(n, mu = 2, sigma = 1, nu = 0.5){
  if (any(mu < 0))  stop(paste("mu must be positive", "\n", ""))  
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", ""))
  if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
  uni<- runif(n = n,0,1)
  r <- qGOLLMax.mixtura(uni,mu =mu, sigma =sigma, nu=nu)
  r
}
#-----------------------------------------------------------------  


dauxiGOLLMax.mixtura <- function(t,mu,sigma,nu,tau){ 
  g <-  dgengamma.orig(t,shape=2,scale=mu,k=3/2)
  G <- pgengamma.orig(t,shape=2,scale=mu,k=3/2)
  fy <- (1-tau)*(nu*sigma*g*(G^(nu*sigma-1))*(1-G^sigma)^(nu-1))/(G^(nu*sigma)+(1-G^sigma)^nu)^2
  fy}



loglog <- function(){
  linkfun<-function(mu) { -log(-log(mu)) }
  linkinv<-function(eta) {
    thresh<-log(-log(.Machine$double.eps))
    eta<-pmin(thresh, pmax(eta, -thresh))
    exp(-exp(-eta)) }
  mu.eta<-function(eta) pmax(exp(-exp(-eta)-eta),.Machine$double.eps)
  valideta<-function(eta) TRUE
  link<-"loglog"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta,
                 name = link), class = "link-gamlss")}
