#####This is R program
##running Poisson INAR function to get estimate and data plug
### Example: INAR(y~x1+x2)
##Nnewton's method is used, \rho's initial value is 0.2
####the model may be not convex, with different initial value, it may return to different optimal result
####To set up your own inital value, run INAR(y~x1+x2,a.init=0.5)

library("ggplot2")

npmtllf<-function(p,y=y,x=matrix(0,nrow=length(y))){ # initialize the parameters
  lp=NULL
  if (length(p)==1){ 
    a<-p[1]
    d<-as.matrix(0)
  } else{
    a<-p[1]
    d<-as.matrix(p[2:length(p)]) 
  }
  y<-as.ts(y)
  x<-as.ts(x)
  
  for(t in 1:length(y)){
    if(t==1) lp[t]=0 else{
      lp[t]=(y[t]-a*y[t-1]-exp(x[t,]%*%d)+a*exp(x[t-1,]%*%d))^2
    }
  }
  lp=-sum(lp)
  lp
}

INAR<-function(formula, a.init=0.2, init.param=NULL, ...)
{ # Set up the model
  u=NULL
  mu=NULL
  model <- model.frame(formula)
  
  # Setup of x and y for the model depends on the number of
  # covariates and the inclusion of the intercept
  
  # Get y and x
  y <- model.extract(model, "response")
  x <- as.matrix(model.matrix(formula))
  xnames <- dimnames(x)[[2]]
  k <- ncol(x)
  
  data_observed=data.frame(time=1:length(y),y=y)
  data_observed$group="Observed"
  
  # Get starting values if needed
  if(is.null(init.param)==TRUE & k>0)
  {
    init.param <- glm(formula, family=poisson())$coefficients
  }
  
  # Now build the correct bounds for the optimizer
  
  bound <- 1000
  if(k>0)
  { lower.bound <- c(0,rep(-bound,k))
  upper.bound <- c(1.1,rep(bound,k))
  } else {
    lower.bound <- c(0)
    upper.bound <- c(1.1)
  }
  
  # Fit the model
  if(k>0)
  {
    fitted.param<-optim(c(a.init,init.param), npmtllf, method = "BFGS",
                        #                            lower=lower.bound, upper=upper.bound,
                        hessian=TRUE,
                        control=list(trace=2, REPORT=5, fnscale=-(length(y))),
                        y=as.matrix(y),x=as.matrix(x))
  } else {
    fitted.param<-optim(c(a.init), npmtllf, method = "L-BFGS-B",
                        lower=lower.bound, upper=upper.bound,
                        hessian=TRUE,
                        control=list(trace=2, REPORT=5, fnscale=-(length(y))),
                        y=as.matrix(y))
  }
  
  covar <- -solve(fitted.param$hessian)
  se<-sqrt(diag(covar))
  # Retrieve the parameters and the llf value
  param<-fitted.param$par
  llf <- fitted.param$value
  aic=-2*llf+2*k
  
  for(t in 1:length(y)){
    if(t==1) u[t]=exp(x[t,]%*%param[2:(k+1)]) else{
      u[t]=param[1]*y[t-1]+exp(x[t,]%*%param[2:(k+1)])-param[1]*exp(x[t-1,]%*%param[2:(k+1)])
    }
  }
  for(t in 1:length(y)){
     mu[t]=exp(x[t,]%*%param[2:(k+1)]) 
  }
  
  data_expected=data.frame(time=1:length(y),y=u)
  data_expected$group="Expected"
  data_plot=rbind(data_observed,data_expected)
  a=ggplot(data_plot,aes(time,y,colour=group))+geom_line(size=1,aes(linetype=group))+
    theme_bw()+xlab("Time")+ylab("Y")+scale_colour_hue("")+scale_linetype_discrete("")+ylim(0,max(data_plot$y)+10)+
    theme(axis.text = element_text(size = 15),
          legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6),
          legend.title=element_text(size=25),legend.text=element_text(size=15),legend.key.size =  unit(0.3, "in"),
          axis.title.x=element_text(size = 15)
    )+xlab("Time")+ylab("")
  rs=y-u
  ls=sum((y-u)^2)
  a=a+annotate("text",x=round((max(data_plot$time)-min(data_plot$time))/5), y= max(data_plot$y)+7,label=paste("AIC:",sprintf("%.2f",round(aic,2))),colour = "red", size = 5.5)
  a=a+annotate("text",x=round((max(data_plot$time)-min(data_plot$time))/5), y= max(data_plot$y)+5,label=paste("SSE:",sprintf("%.2f",round(ls,2))),colour = "red", size = 5.5)
  mylist2=list("coefs"=param,"graph"=a,"AIC"=aic,"Square Error"=ls,"mu"=mu,"residual"=rs,"expected"=u,"standard error"=se)
  return(mylist2)
}
