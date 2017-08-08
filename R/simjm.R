#' Simulate data from a joint model of multiple continuous time-dependent markers and a time-to-event
#'
#' \code{simjm} simulates data from a Cox proportional hazards or semi-parametric additive model with
#' two time-fixed covariates (Z1 and Z2) and three time-dependent covariates (Yij_1,Yij_2,Yij_3).
#' The user can specify various characteristics of these distributions.
#'
#' @export
#'
#' @param n Size of dataset to be generated
#' @param surv_model Model for the time-to-event. Options are "Cox" for Cox proportional hazard model and "Add" for
#' semi-parametric additive hazards model.
#' @param marker_model Model for the multiple markers. Options are "RE" for the correlated random effects model and "PN" for the
#' product normal model.
#' @param MErr Degree of measurement error in the multiple markers. Options are "Low", "Mod" and "High".
#' @param Miss Degree of missing at-risk measurements. Options are "None", "Low" and "High".
#' @param effects Strength of effects of the markers on the hazard. Options are "Null", "Weak" and "Strong".
#' @param corr Degree of marginal pairwise correlations between the markers. Options are "Low", "Mod" and "High".
#' @details The function can be used to generate data from any of the scenarios considered in the the main simulation settings
#' of Moreno-Betancur et al. (2017). See that reference for details.
#' @return A data.frame as required by \code{\link[survtd]{survtd}}. That is, in the long format, with one row per individual and
#' per visit time at which any of the time-dependent covariates were measured, with the corresponding measurements. The dataset also
#' includes a variable that uniquely identifies observations from the same individual; a variable that indicates the timing
#' of each measurement visit; and the fixed variables (time-to-event, event indicator, time-fixed covariates) which are constant across
#' rows of the same individual. The final columns of the dataset (from \code{Xij_1} onwards) are to recover the true values
#' of the markers as per the data generation model for use with function \code{\link[survtd]{simjm_benchmark}}.
#' Specifically, the variables in the dataset are:
#'
#'
#' \describe{
#'   \item{ID}{Unique identifier of observations from the same individual.}
#'   \item{tt}{Time to event, possibly right-censored.}
#'   \item{event}{Indicator of event, with event=1 if an event occurred at tt and event=0 if the individual is censored.}
#'   \item{Z1}{Time-fixed binary covariate.}
#'   \item{Z2}{Time-fixed continuous covariate.}
#'   \item{tj}{Timing of the measurement visit.}
#'   \item{Yij_1}{Measured value of marker 1 at time tj}
#'   \item{Yij_2}{Measured value of marker 2 at time tj }
#'   \item{Yij_3}{Measured value of marker 3 at time tj}
#'   \item{Xij_1}{True value of marker 1 at time tj}
#'   \item{Xij_2}{True value of marker 2 at time tj}
#'   \item{Xij_3}{True value of marker 3 at time tj}
#'   \item{fixed_1}{Time-fixed part of the linear predictor of the linear mixed model from which Yij_1 is generated.}
#'   \item{tim_1}{Time-dependent part of the linear predictor of the linear mixed model from which Yij_1 is generated,
#'                excluding terms for other markers in the case of product-normal model.}
#'   \item{fixed_2}{Time-fixed part of the linear predictor of the linear mixed model from which Yij_2 is generated.}
#'   \item{tim_2}{Time-dependent part of the linear predictor of the linear mixed model from which Yij_2 is generated,
#'                excluding terms for other markers in the case of product-normal model.}
#'   \item{fixed_3}{Time-fixed part of the linear predictor of the linear mixed model from which Yij_3 is generated.}
#'   \item{tim_3}{Time-dependent part of the linear predictor of the linear mixed model from which Yij_3 is generated,
#'                excluding terms for other markers in the case of product-normal model.}
#' }
#'@references
#'Moreno-Betancur M, Carlin JB, Brilleman SL, Tanamas S, Peeters A, Wolfe R (2017). Survival analysis
#'with time-dependent covariates subject to missing data or measurement error: Multiple Imputation for Joint Modeling (MIJM).
#'\emph{Submitted}.
#'
#'@examples
#'
#'   dat<-simjm(n=200,surv_model="Cox",marker_model="PN",
#'              MErr="High",Miss="None",effects="Weak",corr="Low")
#'
#'   head(dat)
#'
#'   dat<-simjm(n=200,surv_model="Add",marker_model="RE",
#'             MErr="Low",Miss="High",effects="Strong",corr="Mod")
#'
#'   head(dat)
#'
#' @import lme4 survival timereg mice
#' @importFrom boot inv.logit
#' @importFrom MASS mvrnorm
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot


simjm<-function(n=200,surv_model="Cox",marker_model="RE",MErr="High",Miss="Low",effects="Weak",
                   corr="Low")
{

######################################## LOAD PARAMETERS ############################################

  #Load all parameters from main set of simulations of Moreno-Betancur et al. based on options above

# TIME-FIXED COVARIATES
pZ1<-0.45
meanZ2<-44
sdZ2<-8.5

# TIME-DEPENDENT MARKERS

# Measurement schedule
J<-	10
Step<-1
tj<-rep(seq(0,J-1,Step),n)
Jt<-length(seq(0,J-1,Step))

# Generation model
alphaZ1_1<-alphaZ1_2<-alphaZ1_3<-(-1.53)
alphaZ2_1<-alphaZ2_2<-alphaZ2_3<-  0.97
alphat_1<-alphat_2<-alphat_3<-	2.5
sdrandI_1<-sdrandI_2<-sdrandI_3<-20
sdrandS_1<-sdrandS_2<-sdrandS_3<-3
if(marker_model=="RE")
{
  alpha0_1<-alpha0_2<-alpha0_3<-	90.13
  rho<-ifelse(corr=="Low",0,ifelse(corr=="Mod",0.5,0.75))
  alphaY2Y1<-alphaY3Y1<-alphaY3Y2<-0
}else
{
  alpha0_1<-90.13
  alpha0_2<-ifelse(corr=="Low",90.13,ifelse(corr=="Mod",75,20))
  alpha0_3<-ifelse(corr=="Low",90.13,ifelse(corr=="Mod",65,10))
  alphaY2Y1<-alphaY3Y1<-ifelse(corr=="Low",0,ifelse(corr=="Mod",0.1,0.5))
  alphaY3Y2<-ifelse(corr=="Low",0,ifelse(corr=="Mod",0.05,0.05))
}
sderror_1<-sderror_2<-sderror_3<-ifelse(MErr=="Low",10,ifelse(MErr=="Mod",20,40))


# TIME-TO-EVENT OUTCOME
v<-2 # (Weibull shape paramter - cannot be changed)

gammaY_1<-ifelse(effects=="Null",0,ifelse(effects=="Weak",0.5,1))

if(surv_model=="Add")
{
  gammaZ1<-	0.01
  gammaZ2<--0.008
  lambda<-ifelse(effects=="Null",0.05,ifelse(effects=="Weak",0.05,0.01))
  gammaY_1<-gammaY_2<-gammaY_3<-gammaY_1/1000

}else
{
  gammaZ1<-	0.6
  gammaZ2<-	0.08
  lambda<-ifelse(effects=="Null",0.0005,ifelse(effects=="Weak",0.00005,0.000007))
  gammaY_1<-gammaY_2<-gammaY_3<-gammaY_1/100

}

# CENSORING
q1<-0.75
q2<-0.85

# MISSING AT-RISK MEASUREMENTS

deltaVis_1<-deltaVis_2<-deltaVis_3<-0.0877
deltaT_1<-deltaT_2<-deltaT_3<-(-0.5)
deltaD_1<-deltaD_2<-deltaD_3<-0.1609
deltaZ1_1<-deltaZ1_2<-deltaZ1_3<-(-0.1113)
deltaZ2_1<-deltaZ2_2<-deltaZ2_3<- 0.0255
deltaYother1_1<-deltaYother1_2<-deltaYother1_3<-0.005
deltaYother2_1<-deltaYother2_2<-deltaYother2_3<-0.005

  if(surv_model=="Add")
{
  if(Miss=="High") delta0_1<-delta0_2<-delta0_3<-ifelse(effects=="Null",1.3,ifelse(effects=="Weak",0,1.3))
  if(Miss=="Low") delta0_1<-delta0_2<-delta0_3<-ifelse(effects=="Null",-0.2,ifelse(effects=="Weak",-1,-0.5))
}else
{
  if(Miss=="High") delta0_1<-delta0_2<-delta0_3<-ifelse(effects=="Null",0.8,ifelse(effects=="Weak",1.3,1.3))
  if(Miss=="Low") delta0_1<-delta0_2<-delta0_3<-ifelse(effects=="Null",-0.7,ifelse(effects=="Weak",-0.5,-0.5))
}

######################################## GENERATE DATA ############################################

ID<-rep(1:n,each=Jt)

# TIME-FIXED COVARIATES
Z1U<-rbinom(n,1,pZ1)
Z1<-rep(Z1U,each=Jt)
Z2U<-rnorm(n,meanZ2,sdZ2)
Z2<-rep(Z2U,each=Jt)

# TIME-DEPENDENT MARKERS
if(marker_model=="RE")
{
  # correlated random intercepts
  dd_I<-sdrandI_1*mvrnorm(n = n, mu=c(0,0,0),
                          Sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),byrow=T,nrow=3,ncol=3))
  a0iU_1<-dd_I[,1]
  a0iU_2<-dd_I[,2]
  a0iU_3<-dd_I[,3]
  # correlated random slopes
  dd_S<-sdrandS_1*mvrnorm(n = n, mu=c(0,0,0),
                          Sigma=matrix(c(1,rho,rho,rho,1,rho,rho,rho,1),byrow=T,nrow=3,ncol=3))
  a1iU_1<-dd_S[,1]
  a1iU_2<-dd_S[,2]
  a1iU_3<-dd_S[,3]
}else
{
  # independent random intercepts and slopes
  a0iU_1<-rnorm(n,0,sdrandI_1)
  a0iU_2<-rnorm(n,0,sdrandI_2)
  a0iU_3<-rnorm(n,0,sdrandI_3)
  a1iU_1<-rnorm(n,0,sdrandS_1)
  a1iU_2<-rnorm(n,0,sdrandS_2)
  a1iU_3<-rnorm(n,0,sdrandS_3)
}

# Marker 1
a0i_1<-rep(a0iU_1,each=Jt)
a1i_1<-rep(a1iU_1,each=Jt)
eij_1<-rnorm(n*Jt,0,sderror_1)
Xij_1<-alpha0_1+alphaZ1_1*Z1+alphaZ2_1*Z2+alphat_1*tj+a0i_1+a1i_1*tj
Yij_1<-Xij_1+eij_1
# Marker 2
a0i_2<-rep(a0iU_2,each=Jt)
a1i_2<-rep(a1iU_2,each=Jt)
eij_2<-rnorm(n*Jt,0,sderror_2)
Xij_2<-alpha0_2+alphaZ1_2*Z1+alphaZ2_2*Z2+alphat_2*tj+a0i_2+a1i_2*tj +alphaY2Y1*Xij_1
Yij_2<-Xij_2+eij_2
# Marker 3
a0i_3<-rep(a0iU_3,each=Jt)
a1i_3<-rep(a1iU_3,each=Jt)
eij_3<-rnorm(n*Jt,0,sderror_3)
Xij_3<-alpha0_3+alphaZ1_3*Z1+alphaZ2_3*Z2+alphat_3*tj+a0i_3+a1i_3*tj +alphaY3Y2*Xij_2+alphaY3Y1*Xij_1
Yij_3<-Xij_3+eij_3
#Check: mean(Yij_1); mean(Yij_2); mean(Yij_3); cor(Yij_1,Yij_2); cor(Yij_1,Yij_3) ; cor(Yij_2,Yij_3)

# TIME-TO-EVENT OUTCOME
if(surv_model=="Add")
{
  ei<-rexp(n)
  K11<-(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1)
  K12<-(alpha0_2+alphaZ1_2*Z1U+alphaZ2_2*Z2U+a0iU_2)+alphaY2Y1*(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1)
  K13<-(alpha0_3+alphaZ1_3*Z1U+alphaZ2_3*Z2U+a0iU_3)+(alphaY3Y1+alphaY3Y2*alphaY2Y1)*(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1)
        +alphaY3Y2*(alpha0_2+alphaZ1_2*Z1U+alphaZ2_2*Z2U+a0iU_2)
  K21<-(alphat_1+a1iU_1)/2
  K22<-(alphat_2+a1iU_2+alphaY2Y1*(alphat_1+a1iU_1))/2
  K23<-(alphat_3+a1iU_3+(alphaY3Y1+alphaY3Y2*alphaY2Y1)*(alphat_1+a1iU_1)+alphaY3Y2*(alphat_2+a1iU_2))/2

  K1<-(gammaY_1*K11+gammaY_2*K12+gammaY_3*K13+gammaZ1*Z1U+gammaZ2*Z2U)
  K2<-(lambda+gammaY_1*K21+gammaY_2*K22+gammaY_3*K23)
  da<-cbind(ei,K1,K2)
  fr<-function(da){
    if(cumhAdd(500,ei=da[1],K1=da[2],K2=da[3])<0) return(202) else
      return(uniroot(cumhAdd,ei=da[1],K1=da[2],K2=da[3],interval=c(0,500))$root)}
}else{
  u<-runif(n)
  K11<-(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1)
  K12<-(alpha0_2+alphaZ1_2*Z1U+alphaZ2_2*Z2U+a0iU_2)+alphaY2Y1*(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1)
  K13<-(alpha0_3+alphaZ1_3*Z1U+alphaZ2_3*Z2U+a0iU_3)+(alphaY3Y1+alphaY3Y2*alphaY2Y1)*(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1)+
    alphaY3Y2*(alpha0_2+alphaZ1_2*Z1U+alphaZ2_2*Z2U+a0iU_2)
  K21<-(alphat_1+a1iU_1)
  K22<-(alphat_2+a1iU_2)+alphaY2Y1*(alphat_1+a1iU_1)
  K23<-(alphat_3+a1iU_3)+(alphaY3Y1+alphaY3Y2*alphaY2Y1)*(alphat_1+a1iU_1)+alphaY3Y2*(alphat_2+a1iU_2)

  K1<-lambda*v*exp(gammaY_1*K11+gammaY_2*K12+gammaY_3*K13+gammaZ1*Z1U+gammaZ2*Z2U)
  K2<-(gammaY_1*K21+gammaY_2*K22+gammaY_3*K23)
  da<-cbind(u,K1,K2)
  fr<-function(da){
    if(cumhCox(500,u=da[1],K1=da[2],K2=da[3])<0) return(202) else
      return(uniroot(cumhCox,u=da[1],K1=da[2],K2=da[3],interval=c(0,500))$root)}
}
ot<-apply(da,1,fr)
minCens<-quantile(ot,prob=q1)
maxCens<-quantile(ot,prob=q2)
cens<-runif(n,minCens,maxCens)
tt<-ot*(ot<=cens)+cens*(ot>cens)
event<-1*(ot<=cens)+0*(ot>cens)
# Check: sum(ot>cens)/n; hist(ot); hist(tt)

# DATASET WITH TIME-DEPENDENT MARKERS (dat): one row per individual per measurement
dat<-as.data.frame(cbind(ID=ID,tt=rep(tt,each=Jt),event=rep(event,each=Jt),Z1,Z2,
                         fixed_1=rep(alpha0_1+alphaZ1_1*Z1U+alphaZ2_1*Z2U+a0iU_1,each=Jt),
                         tim_1=rep(alphat_1+a1iU_1,each=Jt),
                         fixed_2=rep(alpha0_2+alphaZ1_2*Z1U+alphaZ2_2*Z2U+a0iU_2,each=Jt),
                         tim_2=rep(alphat_2+a1iU_2,each=Jt),
                         fixed_3=rep(alpha0_3+alphaZ1_3*Z1U+alphaZ2_3*Z2U+a0iU_3,each=Jt),
                         tim_3=rep(alphat_3+a1iU_3,each=Jt),
                         Xij_1,Yij_1,Xij_2,Yij_2,Xij_3,Yij_3,tj))
dat<-dat[dat$tt>=dat$tj,]

# GENERATE MISSING AT RISK MEASUREMENTS
flag<-0
if(Miss!="None")
{
  # Marker 1
  dat$Proba_1<-(dat$tj!=0)*inv.logit(delta0_1+deltaT_1*dat$tt+deltaD_1*dat$event+
                                       deltaZ1_1*dat$Z1+deltaZ2_1*dat$Z2+deltaVis_1*dat$tj+
                                       deltaYother1_1*dat$Xij_2+deltaYother2_1*dat$Xij_3)
  flag<-0
  while(flag==0)
  {
    dat$mis_1<-rbinom(nrow(dat),1,dat$Proba_1)
    gm<-sum(dat$mis_1)/nrow(dat)

    if(Miss=="Low")
    {if(gm>=0.15&gm<=0.25) flag<-1 }else
      if(Miss=="High")
      {if(round(gm,2)>=0.35&round(gm,1)<=0.5) flag<-1 }
  }
  dat[dat$mis_1==1,"Yij_1"]<-NA

  # Marker 2
  dat$Proba_2<-(dat$tj!=0)*inv.logit(delta0_2+deltaT_2*dat$tt+deltaD_2*dat$event+
                                       deltaZ1_2*dat$Z1+deltaZ2_2*dat$Z2+deltaVis_2*dat$tj+
                                       deltaYother1_2*dat$Xij_1+deltaYother2_2*dat$Xij_3)

  flag<-0
  while(flag==0)
  {
    dat$mis_2<-rbinom(nrow(dat),1,dat$Proba_2)
    gm<-sum(dat$mis_2)/nrow(dat)

    if(Miss=="Low")
    {if(gm>=0.15&gm<=0.25) flag<-1 }else
      if(Miss=="High")
      {if(round(gm,2)>=0.35&round(gm,1)<=0.5) flag<-1 }
  }

  dat[dat$mis_2==1,"Yij_2"]<-NA


  dat$Proba_3<-(dat$tj!=0)*inv.logit(delta0_3+deltaT_3*dat$tt+deltaD_3*dat$event+
                                       deltaZ1_3*dat$Z1+deltaZ2_3*dat$Z2+deltaVis_3*dat$tj+
                                       deltaYother1_3*dat$Xij_1+deltaYother2_3*dat$Xij_2)

  flag<-0
  while(flag==0)
  {
    dat$mis_3<-rbinom(nrow(dat),1,dat$Proba_3)
    gm<-sum(dat$mis_3)/nrow(dat)

    if(Miss=="Low")
    {if(gm>=0.15&gm<=0.25) flag<-1 }else
      if(Miss=="High")
      {if(round(gm,2)>=0.35&round(gm,1)<=0.5) flag<-1 }
  }

  dat[dat$mis_3==1,"Yij_3"]<-NA
}
dat<-dat[,c("ID","tt","event","Z1","Z2","tj","Yij_1", "Yij_2","Yij_3", "Xij_1",
            "Xij_2","Xij_3","fixed_1","tim_1", "fixed_2","tim_2","fixed_3","tim_3")]

return(dat)
}








