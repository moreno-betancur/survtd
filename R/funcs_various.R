#' Cumulative hazard for Cox model with Weibull baseline.
#'
#' Internal function used to generate data from Cox model.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

cumhCox<-function(tt,u,K1,K2){
  if(K2!=0)return(K1*((tt*exp(K2*tt)/K2)-(exp(K2*tt)/(K2^2))+(1/(K2^2)))+log(u)) else
    return(K1*(tt^2)/2+log(u))
}

#' Cumulative hazard for additive model with Weibull baseline.
#'
#' Internal function used to generate data from additive model.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

cumhAdd<-function(tt,ei,K1,K2){
  return(K1*tt+K2*(tt^2)-ei)
}

#' Identifies time-dependent covariates in the model
#'
#' Specifies which of the regressors in the model are time-dependent.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

td<-function(x)x

#' Derive modified event indicator for imputation
#'
#' Used to derive the modified event indicator on the dataset in long format.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot
#'
evMI<-function(dataset,visit.time,ename)
{ dataset$event2<-0
dataset[which(dataset[,visit.time]==max(dataset[,visit.time])),"event2"]<- dataset[which(dataset[,visit.time]==max(dataset[,visit.time])),ename]
return(dataset)
}

#' Apply last observation carried forward to long-format dataset
#'
#' In the dataset in long format, uses the last observation carried forward (LOCF) approach to fill in missing values.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

locf<-function(dataset,marker)
{ if(nrow(dataset)>1)
{for(k in 2:(nrow(dataset)))
{ if(is.na(dataset[k,marker]))  dataset[k,marker]<-dataset[k-1,marker]
}}
  return(dataset)
}

#' Summary method for semi-parametric additive model fitted with aalen function
#'
#' Extracts hazard differences, confidence intervals and p-values for model fitted using \code{\link[timereg]{aalen}} function.
#' @export
#' @keywords internal
#' @importFrom stringr str_replace
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

summ.aalen<-function(fit)
{
  Estimate<-t(fit$gamma)
  SE<-sqrt(diag(fit$var.gamma))
  CILow<-Estimate-1.96*SE
  CIUpp<-Estimate+1.96*SE
  Pvalue<-2*(1-pnorm(abs(Estimate/SE)))
  out<-data.frame(t(Estimate),SE,t(CILow),t(CIUpp),t(Pvalue))
  names(out)<-c("HD","SE","CIlow","CIupp","p-value")
  nam<-row.names(out)
  nam<-str_replace(nam,"const\\(","")
  nam<-str_replace(nam,"\\)","")
  row.names(out)<-nam
  return(out)
}

#' Summary method for Cox model fitted with coxph function
#'
#' Extracts hazard ratios, confidence intervals and p-values for model fitted using \code{\link[survival]{coxph}} function.
#' @export
#' @keywords internal
#' @import survival
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

summ.cox<-function(fit)
{
  out<-data.frame(summary(fit)$coef[,c("coef")],summary(fit)$coef[,c("se(coef)")],
                  log(summary(fit)$conf.int[,3:4]),summary(fit)$coef[,5])
  names(out)<-c("logHR","SE","CIlow","CIupp","p-value")
  return(out)
}


#' Pool estimes from multiple imputed datasets
#'
#' Combines the point and variance estimates obtained from each imputed dataset into the final multiple imputation estimate using Rubin's rules.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

combMI<-function(coefmi,varmi,M)
{est<-mean(coefmi)
varinter<-(1/(M-1))*sum((coefmi-est)^2)
varintra<-mean(varmi)
se<-sqrt(varintra+(1+1/M)*varinter)
r<-(1+1/M)*varinter/varintra
vvv<-(M-1)*(1+1/r)^2
tal<-qt(0.025,df=vvv,lower.tail=F)
pval<-2*(1-pt(q=abs(est/se),df=vvv))
return(c(est,se,est-tal*se,est+tal*se,pval))}



#' Summarise results from MICEa and MICEb methods
#'
#' Extracts final multiple imputation estimates of hazard ratios or differences, according to the model, and confidence intervals and p-values
#' for models fitted using MICEa or the MICEb approach.
#' @export
#' @keywords internal
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

summ.survtd<-function(coefs,vars,model,nam,M)
{
  out<-data.frame()
  for(j in 1:ncol(coefs))
    out<-rbind(out,combMI(coefs[,j],vars[,j],M))
  if(model=="Add")
    names(out)<-c("HD","SE","CIlow","CIupp","p-value")
  else
  {
    names(out)<-c("logHR","SE","CIlow","CIupp","p-value")
  }
  row.names(out)<-nam
  return(out)
}



