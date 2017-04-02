#' Fit survival hazard models with time-dependent covariates
#'
#' \code{survtd} fits Cox proportional hazards and semi-parametric additive hazards models
#' with time-fixed covariates (any type) and time-dependent covariates either with a multiple imputation-based two-stage
#' joint modelling approach to deal with missing data and measurement error (only continuous time-dependent covariates),
#' a naive two-stage approach (only continuous time-dependent covariates) or the LOCF approach (any type).
#'
#'
#' @export
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on the right as regressors.
#' The response must be a survival object as returned by the Surv function of type "right" (other types are not
#' supported at this stage). Time-dependent regressors are specified by the wrapper \code{td()} (see \strong{Examples} section below).
#' This development version does not support interactions nor terms constructed using I(). The user thus needs to create the
#' necessary interaction/derived variables in the dataset before using this function.
#' @param data A data.frame in which to interpret the variables named in the formula. The dataset needs to be in long format, with
#' one row per individual and per visit time at which any of the time-dependent covariates were measured, with the corresponding measurements.
#' The dataset must also include a variable that uniquely identifies observations from the same individual (see parameter \code{id} below);
#' a variable that indicates the timing of each measurement visit (see parameter \code{visit.time} below); and other fixed variables
#' (time-to-event, event indicator, time-fixed covariates) which are constant across rows of the same individual.
#' @param id The name of the variable in \code{data} that uniquely identifies observations from the same individual.
#' @param visit.time The name of the variable in \code{data} that indicates the timing of each visit.
#' @param model Indicates which hazard model to fit. Options are "Cox" for the Cox proportional hazard model and "Add" for the semi-parametric
#' additive hazards model.
#' @param method Indicates which method to use to fit the model. Options are "MICEa", "MICEb", "LOCF" and "Naive2S". Both "MICEa" and "MICEb" use a two-stage
#' joint modelling approach based on multiple imputation to incorporate time-varying continuous covariates (see \strong{Details} section below).
#' Method "LOCF" uses the Last Observation Carried Forward (LOCF) approach. Method "Naive2S" is the naive two-stage approach (see \strong{Details} section below).
#' @param M Number of imputations to perform for methods "MICEa" and "MICEb" (ignored if method="LOCF" or "Naive2S").
#' @param G Number of iterations to perform in multiple imputation by chained equations (MICE) algorithm for methods "MICEa" and "MICEb" (ignored if method="LOCF" or "Naive2S").
#' @param time.trend Formula object with empty left hand side and right hand side expressing a polynomial of "x" which determines the
#' way time is modeled in the fixed effects part of the linear mixed model for each time-depdent marker (ignored if method="LOCF"). Default is a linear trend
#' (i.e. include \code{visit.time} as predictor in the model). The random effects part includes a random intercept and slope.
#' Future development plans are to allow a different \code{time.trend} argument for each marker, the possibility to include
#' natural cubic splines in this argument and more general random effects structures.
#'
#' @details The \code{survtd} function can be used to fit the Cox proportional hazards or
#'   semi-parametric additive hazard models with time-depedent covariates.
#'
#'   Methods "MICEa" and "MICEb" can be used to fit models with time-fixed
#'   covariates of any type and multiple continuous time-dependent covariates.
#'   To deal with the discrete-time observation of the time-varying covariates,
#'   particularly measurement error and missing data, a two-stage joint modeling approach
#'   is used that is based on multiple imputation. Details are provided in Moreno-Betancur
#'   et al. (2016), but briefly, in Stage 1 the true (error-corrected) values of the
#'   time-dependent covariates at each event time at which an individual is at risk are multiply
#'   imputed by drawing iteratively from smoothed trajectories based on interdependent
#'   linear mixed models using Multiple Imputation by Chained Equations (MICE).
#'   An adaptation of the \code{mice} function from the \emph{mice} package is used for this step.
#'   In Stage 2, the time-to-event model is fitted to each of the imputed datasets
#'   (see below) and estimates are pooled using Rubin's MI formulas.
#'
#'   The two MICE methods differ in the way information about the event occurrence is
#'   included in the imputation models: the "MICEb" method includes the
#'   the event indicator, while the "MICEa" approach is based on a more refined
#'   approximation including a modified event indicator (see Moreno-Betancur et al. 2016).
#'   Hence, the approach "MICEa" is generally preferable.
#'
#'   The "Naive2S" method can be used to fit models with time-fixed covariates of any type
#'   and multiple continuous time-dependent covariates using the naive two-stage approach.
#'   This method singly imputes each continuous marker at each event time from its estimated
#'   trajectory obtained from a linear mixed model including fixed and random effects for time
#'   and fixed effects for the time-fixed covariates appearing in the formula argument. The
#'   method thus ignores: the uncertainty in these imputed values, the interrelations between
#'   the time-dependent markers and the relation between the time-dependent markers and the
#'   time-to-event process (see Moreno-Betancur et al. 2016).
#'
#'   The "LOCF" method can be used with any type of time-fixed and time-varying covariates
#'   (categorical or continuous). This approach uses the last available measurement of the
#'   time-varying covariate to singly impute its value at each event time at which an
#'   individual is at risk. It can perform very poorly if the observations are not
#'   synchronously updated across individuals (e.g. due to missing data) and if the
#'   time-varying covariates are measured with error (see Moreno-Betancur et al. 2016).
#'
#'   Once the values of the time-dependent covariates are imputed at each of the event-times at which
#'   the individual is at risk according to either of the four methods,
#'   the function uses \code{coxph} from the \emph{survival} package to fit the Cox model and
#'   \code{aalen} from the \emph{timereg} package to fit the additive model.
#'
#'@return This development version returns a data frame with regression coefficient estimates for each
#'covariate based on the method chosen, along with standard errors, 95\%  confidence intervals and p-values.
#'This will change in future versions when a proper class of objects and summary and other such methods are developed.
#'
#'@references Moreno-Betancur M, Carlin JB, Brilleman SL, Tanamas S, Peeters A, Wolfe R (2016). Survival analysis with
#'time-dependent covariates subject to measurement error and missing data: Two-stage joint model using multiple imputation. Submitted.
#'
#'
#'@examples
#'
#'   ## Example with additive model ##
#'
#'   dat<-simjm(n=200,surv_model="Add",marker_model="RE",
#'             MErr="Low",Miss="High",effects="Strong",corr="Mod")
#'
#'   survtd(Surv(time=tt,event=event)~Z1+Z2+td(Yij_1)+td(Yij_2)+td(Yij_3),
#'          data=dat,  id="ID", visit.time="tj", model="Add",
#'          method="MICEa", M=5, G=5,time.trend=as.formula("~x+I(x^2)"))
#'
#'   survtd(Surv(time=tt,event=event)~Z1+Z2+td(Yij_1)+td(Yij_2)+td(Yij_3),
#'          data=dat,  id="ID", visit.time="tj", model="Add",
#'          method="Naive2S", M=5, G=5,time.trend=as.formula("~x+I(x^2)"))
#'
#'   survtd(Surv(time=tt,event=event)~Z1+Z2+td(Yij_1)+td(Yij_2)+td(Yij_3),
#'          data=dat,  id="ID", visit.time="tj", model="Add",
#'          method="LOCF", M=5, G=5)
#'
#'   simjm_benchmark(dat,surv_model="Add",marker_model="RE",corr="Mod")
#'
#' @import lme4 survival timereg mice
#' @importFrom stringr str_replace
#' @importFrom ipw tstartfun
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot model.matrix


survtd<-function(formula, data, id,visit.time,model="Cox",method="MICEa",M=5,G=5, time.trend=as.formula("~x"))
{
##################################### Preliminaries and checks #####################################

  call <- match.call()
  m <- match.call(expand.dots = F)
  if (!model %in% c("Cox", "Add"))
    stop("model must be \"Cox\" or \"Add\"")
  if (!method %in% c("LOCF", "MICEa", "MICEb", "Naive2S"))
    stop("method must be \"LOCF\", \"MICEa\", \"MICEb\" or \"Naive2S\" ")
  special <- c("td")
  Terms <- terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m$id <- m$visit.time <- m$G <- m$M <- m$model <- m$method <-m$time.trend<- NULL
  m$na.action <- NULL
  m <- eval(m, parent.frame())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  S <- model.extract(m, "response")
  if (!inherits(S, "Surv"))
    stop("Response must be a survival object")


# Recover names of all variables in model
  tname <- as.character(Terms[[2]][2])
  ename <- as.character(Terms[[2]][3])
  fnames <- attr(Terms, "term.labels")[-(attr(Terms, "specials")$td -
                                           1)]
  tdnames <- attr(Terms, "term.labels")[(attr(Terms, "specials")$td -
                                           1)]
  tdnames <- str_replace(tdnames, "td\\(", "")
  tdnames <- str_replace(tdnames, "\\)", "")
  if (any(apply(data[,names(data) %in% tdnames], 2, class) != "numeric") & method != "LOCF")
    stop("For methods MICEa, MICEb and Naive2S, time-dependent variables must be of class \"numeric\"")

##################################### Data preparation ##############################################


# Recover dataset with one row per individual and time-fixed variables to get NA cumulative hazard
# estimator and interactions with time-fixed covariates
  datU <- unique(cbind(data[, c(id, tname, ename, fnames)]))
  cumhz <- basehaz(coxph(Surv(datU[, tname], datU[, ename]) ~
                           1))
  time <- datU[, tname]
  cumhz <- merge(as.data.frame(time), cumhz, by = c("time"),
                 all.x = TRUE, sort = F, suffixes = c("", "y"))[, "hazard"]
  datU$cumhz <- cumhz
  intnames<-vector()
  for (Z in fnames) {
    if(is.factor(datU[,Z]))
      mint<-as.data.frame(model.matrix(as.formula(paste("~",Z,"-1")),datU)[,-1])
    else
      mint<-as.data.frame(model.matrix(as.formula(paste("~",Z,"-1")),datU))

    naint<-paste("int", Z, 1:ncol(mint),sep = "")
    names(mint)<-naint
    datU<-cbind(datU,mint*datU$cumhz)
    intnames<-c(intnames,naint)}

# Add these variables to original dataset

  data <- merge(data, datU[, c(id, "cumhz", intnames)], by = id, all.x = TRUE, suffixes = c("x",  ""))
  data <- do.call("rbind", by(data, data[, id], evMI, visit.time,
                              ename, simplify = F))

# Base dataset for time-to-event analyses (datSurv): one row per individual per observed event-time
  times <- sort(unique(c(datU[, tname][datU[, ename] != 0])))
  ID <- rep(datU[, id], each = length(times))
  tt <- rep(datU[, tname], each = length(times))
  fuptime <- rep(times, nrow(datU))
  datSurv <- data.frame(ID, tt, fuptime)
  datSurv <- datSurv[datSurv$fuptime <= datSurv$tt, ]
  datSurv$tstart <- tstartfun(ID, fuptime, datSurv)
  names(datSurv)[1] <- id
  datSurv$tstart[datSurv$tstart == -1] <- 0
  datSurv <- merge(datSurv, datU[, !names(datU) %in% c(tname)],
                   by = id, all.x = TRUE)
  datSurv$event2 <- with(datSurv, ifelse(fuptime != tt, 0,
                                         get(ename)))
  datSurv[, visit.time] <- ifelse(floor(datSurv$fuptime)%%2 ==
                                    0, floor(datSurv$fuptime), floor(datSurv$fuptime) - 1)
  datSurv[, visit.time] <- ifelse(datSurv[, visit.time] > max(data[,
                                                                   visit.time]), max(data[, visit.time]), datSurv[, visit.time])

  datLOCF <- data
  for (mark in tdnames) datLOCF <- do.call("rbind", by(datLOCF,
                                                       datLOCF[, id], locf, marker = mark, simplify = F))
  datSurv <- merge(datSurv, datLOCF[, c(id, visit.time, tdnames)],
                   by = c(id, visit.time), all.x = TRUE, sort = F, suffixes = c("","y"))

######################################## ANALYSIS ######################################################

if(method=="LOCF")     #################  LOCF approach #####################
{
  if (model == "Cox") {
    fit <- coxph(as.formula(paste("Surv(time=tstart,time2=fuptime,event=event2)~",
                                  paste(c(tdnames, fnames), collapse = "+"))),
                 datSurv)
    res <- summ.cox(fit)
  }
  else {
    fit <- aalen(as.formula(paste("Surv(time=tstart,time2=fuptime,event=event2)~",
                                  paste(paste("const(", c(tdnames, fnames), ")",
                                              sep = ""), collapse = "+"))), data = datSurv,
                 robust = F, n.sim = 0)
    res <- summ.aalen(fit)
  }
}else{                ################   Two-stage approaches #################

#### Create dataset for 2-stage approaches: constructed by stacking two datasets

# 1. Dataset with observed marker values from which the imputation model will be built
  datl <- data
  datl$tstart <- (-1)
  datl$fuptime <- (-1)
  datl$ind <- 1
  datl <- datl[, c(id, "tstart", "fuptime", visit.time,
                   "ind", tname, "cumhz", "event2", ename, fnames,  intnames, tdnames)]

# 2. Dataset with one row per individual per observed event-time, with all marker values missing
#    This will enable imputation of the markers at the exact event-times
dat2S <- datSurv
dat2S <- merge(datSurv, datLOCF[, c(id, visit.time, tdnames)],
               by = c(id, visit.time), all.x = TRUE, sort = F, suffixes = c("",
                                                                            "y"))
dat2S[, visit.time] <- dat2S$fuptime
dat2S[, tdnames] <- NA
dat2S$ind <- 2
dat2S[, tname] <- dat2S$tt
dat2S <- dat2S[, c(id, "tstart", "fuptime", visit.time,
                   "ind", tname, "cumhz", "event2", ename, fnames, intnames, tdnames)]


# Stack the two datasets
dat2S<-rbind(datl,dat2S)

#prepare time-trend variables
if(time.trend==as.formula("~x"))
  vnames<-NULL else
  {
    datvis<-as.data.frame(dat2S[,visit.time])
    names(datvis)<-"x"
    matvis<-as.data.frame(model.matrix(time.trend,datvis)[,-1])
    matvis<-as.data.frame(matvis[,!names(matvis)%in%"x"])
    vnames<-paste("time.term",1:ncol(matvis),sep="")
    names(matvis)<-vnames
    dat2S<-cbind(dat2S, matvis)}

if(method=="MICEa"|method=="MICEb")
{          ######### Proposed 2-stage approach based on MICE #####################

# Stage 1: Multiply impute error-free value of markers at all observed event times

# Initialise and then set predictor matrix to feed to mice, where non-zero entries indicate
# that the variable of the column is included in linear mixed imputation model of row variable
# Specifically, the following codes are for assigning values to each entry:
# 0  = indicates that variable does not enter the imputation model
# -2 = indicates subject ID variable
# 1  = indicates that variable has a fixed effect only
# 2  = indicates that variable has fixed effect and a subject-specific random effect
  ini <- mice(dat2S, maxit = 0)
  pred <- ini$pred
  ev <- ifelse(method == "MICEa", "event2", ename)
  ti <- "cumhz"
  for (mark in tdnames) {
    pred[mark, ] <- 0
    pred[mark, id] <- (-2)
    pred[mark, visit.time] <- 2
    pred[mark, c(fnames, intnames,
                 setdiff(tdnames, mark), ev, ti,vnames)] <- 1
  }
  mm <- rep("", dim(pred)[1])
  mm[dimnames(pred)[[1]] %in% tdnames] <- "lmm"

  impMICE <- mice2(dat2S, m = M, method = mm, predictorMatrix = pred,
                   maxit = G, printFlag = F)

# Stage 2: Fit the survival model to each imputed dataset and combine the estimates using Rubin's rules

  coefs <- data.frame()
  vars <- data.frame()
  for (m in 1:M) {
    datcom <- complete(impMICE, m)
    datcom <- datcom[datcom$ind == 2, ]
    if (model == "Cox") {
      fit <- coxph(as.formula(paste("Surv(time=tstart,time2=fuptime,event=event2)~",
                                    paste(c(tdnames, fnames), collapse = "+"))),
                   data = datcom)
      coefs <- rbind(coefs, summary(fit)$coef[, c("coef")])
      vars <- rbind(vars, (summary(fit)$coef[, c("se(coef)")])^2)
      nam <- attr(fit$coef, "names")
    }
    else {
      fit <- aalen(as.formula(paste("Surv(time=tstart,time2=fuptime,event=event2)~",
                                    paste(paste("const(", c(tdnames, fnames),
                                                ")", sep = ""), collapse = "+"))), data = datcom,
                   robust = F, n.sim = 0, silent = 0)
      coefs <- rbind(coefs, t(fit$gamma))
      vars <- rbind(vars, diag(fit$var.gamma))
      nam <- attr(fit$gamma, "dimnames")[[1]]
      nam <- str_replace(nam, "const\\(", "")
      nam <- str_replace(nam, "\\)", "")
    }
  }
  res <- summ.survtd(coefs, vars, model, nam, M)
}else if(method=="Naive2S"){  ############## Naive two-stage approach #############
  datNaive <- dat2S
  # Stage 1: Singly impute error-free value of markers at all observed event times

  for (mark in tdnames) {
    fitNaive <- lmer(as.formula(paste(mark, "~1+",
                                      paste(fnames, collapse = "+"), "+",paste(vnames, collapse = "+"),"+", visit.time,
                                      "+(1+", visit.time, "|", id, ")")), data = datNaive)
    datNaive[, mark] <- lme4:::predict.merMod(fitNaive,
                                              newdata = datNaive, type = "response", re.form = NULL)
  }
  # Stage 2: Fit the survival model to singly imputed dataset

  datNaive <- datNaive[datNaive$ind == 2, ]
  if (model == "Cox") {
    fit <- coxph(as.formula(paste("Surv(time=tstart,time2=fuptime,event=event2)~",
                                  paste(c(tdnames, fnames), collapse = "+"))),
                 data = datNaive)
    res <- summ.cox(fit)
  }else {
    fit <- aalen(as.formula(paste("Surv(time=tstart,time2=fuptime,event=event2)~",
                                  paste(paste("const(", c(tdnames, fnames), ")",
                                              sep = ""), collapse = "+"))), data = datNaive,
                 robust = F, n.sim = 0, silent = 0)
    res <- summ.aalen(fit)
  }

}
}

  return(list(Call=call,Results=res))
}