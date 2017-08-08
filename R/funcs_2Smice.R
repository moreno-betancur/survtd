#' Draw imputations from a linear mixed model
#'
#' Implements the imputation method "lmm" to specify in mice2 function to draw imputations from linear mixed imputation models.
#' @export
#' @keywords internal
#' @import lme4 mice
#' @import Matrix
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

mice.impute.lmm <-
  function (ry, namesX, namesZ, nameG, nam, data)
  {
        ### CREATE MODEL FORMULA ###
        if(length(nameG) == 0)
          stop(paste("No grouping factor to include random effects in model."))

        if (length(namesZ) == 0) {
          formulaGLMER = as.formula(paste(nam,"~1+",paste(namesX,collapse="+"),
                                          paste("+(1|",nameG,")")))
        } else if (length(namesZ) != 0)
          formulaGLMER = as.formula(paste(nam,"~1+",paste(namesX,collapse="+"),
                                          paste("+(1+",paste(namesZ,collapse="+"),"|",nameG,")")))

        ### FIT LMM FROM AVAILABLE CASES ###
        fit <- lme4::lmer(formulaGLMER, data[ry==1,])
        ids <- sort(unique(unlist(data[ry==1, nameG], FALSE, FALSE)))

        ### RETRIEVE FIXED EFFECTS ESTIMATES ###
        Beta    <- lme4::fixef(fit)
        VarBeta <- t(Matrix::chol(as.matrix(vcov(fit))))

        ### RETRIEVE RANDOM EFFECTS PREDICTIONS AND CONDITIONAL VAR-COV MATRIX ###
        Bmat    <- lme4::ranef(fit, condVar = TRUE)[[1]]  # nLevels x nRanef matrix
        nLevels <- nrow(Bmat)
        nRanef  <- ncol(Bmat)
        B       <- unlist(Bmat, FALSE, FALSE)
        VarB    <- attr(Bmat, "postVar", exact = TRUE)

        ### DRAW FIXED AND RANDOM EFFECTS FROM THEIR ASYMPTOTIC DISTRIBUTION ###
        beta.star <- Beta + VarBeta %*% rnorm(ncol(VarBeta))
        B.disturbance <- as.vector(t(
          do.call(cbind, lapply(seq(dim(VarB)[3]),
                                function(x) {
                                  L <- t(Matrix::chol(VarB[ , , x]))
                                  B.disturbance <- L %*% rnorm(ncol(L))
                                }))))
        b.star    <- B + B.disturbance

        ### OBTAIN MODEL MATRICES FOR IMPUTATION ###
        data2 <- data
        data2[is.na(data2[, nam]), nam] <- 0
        ids2  <- sort(unique(unlist(data2[!ry, nameG], FALSE, FALSE)))
        if (!identical(ids, ids2))
          stop(paste0("Trying to impute for individuals who did not have ",
                      "random effect estimates from available case data."))

        formulaGLMER2 = as.formula(paste(nam,"~1+",paste(namesX,collapse="+")))
        X <- model.matrix(formulaGLMER2, data2[!ry, ])

        formulaGLMER3 = as.formula(paste(nam,"~-1+as.factor(",nameG,")+as.factor(",nameG,"):",paste(namesZ,collapse="+")))
        Z <- Matrix::sparse.model.matrix(formulaGLMER3, data2[!ry, ])

        ### VALUES TO IMPUTE MISSING OUTCOMES = LINEAR PREDICTOR  ###
        fixedpart <- X %*% beta.star
        vec <- fixedpart + Z %*% b.star
        return(vec@x)
  }
environment(mice.impute.lmm)<-environment(mice)

#' Internal function for multiple imputation by chained equations allowing for "lmm" imputation method.
#'
#' Modified version of \code{sampler} function from the \code{mice} package to allow for "lmm" imputation method.
#' @export
#' @keywords internal
#' @import mice
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot


sampler2 <- function(p, data, m, imp, r, visitSequence, fromto, printFlag, ...)
  # The sampler controls the actual Gibbs sampling iteration scheme This function is called by mice and mice.mids
  #
  # Authors: S van Buuren, K Groothuis-Oudshoorn Copyright (c) 1999-2008 TNO Quality of Life
{
  ## set up array for convergence checking
  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  if (maxit > 0)
    chainVar <- chainMean <- array(0, dim = c(length(visitSequence), maxit, m), dimnames = list(dimnames(data)[[2]][visitSequence],
                                                                                                1:maxit, paste("Chain", 1:m))) else chainVar <- chainMean <- NULL

  ## THE ITERATION MAIN LOOP: GIBBS SAMPLER
  if (maxit < 1)
    iteration <- 0 else {
      if (printFlag)
        cat("\n iter imp variable")
      for (k in from:to) {
        #begin k loop : iteration loop
        iteration <- k
        for (i in 1:m) {
          #begin i loop    : repeated imputation loop
          if (printFlag)
            cat("\n ", iteration, " ", i)

          ## fill the data with the last set of imputations
          for (j in visitSequence) p$data[!r[, j], j] <- imp[[j]][, i]

          ## augment the data with the actual dummy variables
          for (j in setdiff(p$visitSequence, visitSequence)) {
            cat.columns <- p$data[, p$categories[j, 4]]
            p$data[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- matrix((model.matrix(~cat.columns - 1)[, -1]),
                                                                                  ncol = p$categories[p$categories[j, 4], 2], nrow = nrow(p$data))
          }

          ## iterate once over the variables of the augmented model

          for (j in p$visitSequence) {
            theMethod <- p$method[j]
            vname <- dimnames(p$data)[[2]][j]

            ## store current state
            oldstate <- get("state", pos = parent.frame())
            newstate <- list(it = k, im = i, co = j, dep = vname, meth = theMethod, log = oldstate$log)
            assign("state", newstate, pos = parent.frame(), inherits = TRUE)

            if (printFlag & theMethod != "dummy")
              cat(" ", vname)
            if (theMethod != "" & (!is.passive(theMethod)) & theMethod != "dummy") {
              # for a true imputation method
              if (substring(tolower(theMethod), 1, 2) != "2l" &
                    substr(theMethod,1,3) != "lmm") {
                # RJ: for an non-multilevel imputation method
                # RB: formula-based  specification
                if (! is.null(p$form) && nchar(p$form[j])>0) {
                  myform <- paste(p$form[j], "0", sep="+")
                  x <- model.matrix(formula(myform), p$data)
                } else
                  x <- p$data[, p$predictorMatrix[j, ] == 1, drop = FALSE]
                y <- p$data[, j]
                ry <- r[, j]
                nam <- vname
                if (k == 1)
                  check.df(x, y, ry, ...)  # added 31/10/2012, throw warning for n(obs) < p case
                f <- paste("mice.impute", theMethod, sep = ".")
                keep <- remove.lindep(x, y, ry, ...)
                x <- x[, keep, drop = FALSE]
                imp[[j]][, i] <- do.call(f, args = list(y, ry, x, ...))
              } else if (substring(theMethod, 1, 2) == "2l"){
                # for a multilevel imputation method
                predictors <- p$predictorMatrix[j, ] != 0
                # RB: formula-based specification
                if (! is.null(p$form) && nchar(p$form[j])>0) {
                  myform <- paste(p$form[j], "0", sep="+")
                  x <- model.matrix(formula(myform), p$data)
                } else
                  x <- p$data[, predictors, drop = FALSE]
                y <- p$data[, j]
                ry <- r[, j]
                type <- p$predictorMatrix[j, predictors]
                nam <- vname
                if (k == 1)
                  check.df(x, y, ry, ...)  # added 31/10/2012, throw warning for n(obs) < p case
                f <- paste("mice.impute", tolower(theMethod), sep = ".")
                keep <- remove.lindep(x, y, ry, ...)
                x <- x[, keep, drop = FALSE]
                type <- type[keep]
                imp[[j]][, i] <- do.call(f, args = list(y, ry, x, type, ...))
              } else if(substring(theMethod, 1, 3) == "lmm"){
                predictors_fixed <- p$predictorMatrix[j, ]>0
                predictors_random<- p$predictorMatrix[j, ] ==2
                group <- p$predictorMatrix[j, ] ==-2
                x <- p$data[, predictors_fixed, drop = FALSE]
                z <- p$data[, predictors_random, drop = FALSE]
                y <- p$data[, j]
                ry <- r[, j]
                nam <- dimnames(p$data)[[2]][j]
                namesX1<-dimnames(p$data)[[2]][predictors_fixed]
                namesZ1<-dimnames(p$data)[[2]][predictors_random]
                nameG<-dimnames(p$data)[[2]][group]
                nameG<-strsplit(as.character(nameG[1]),split="\\.")[[1]][1]
                f <- paste("mice.impute", theMethod, sep = ".")
                keep <- remove.lindep(x, y, ry)
                x<- x[, keep, drop = FALSE]
                namesX<-namesX1[keep]
                namesZ<-intersect(namesX,namesZ1)
                z <- z[, namesZ, drop = FALSE]
                imp[[j]][, i] <- mice.impute.lmm(ry, namesX, namesZ, nameG, nam,p$data)
              }
              p$data[!r[, j], j] <- imp[[j]][, i]
            } else if (is.passive(theMethod)) {
              imp[[j]][, i] <- model.frame(as.formula(theMethod), p$data[!r[, j], ])  #RJ - FIXED passive imputation: as.formula()
              p$data[!r[, j], j] <- imp[[j]][, i]
            } else if (theMethod == "dummy") {
              ## FEH
              cat.columns <- p$data[, p$categories[j, 4]]
              p$data[, (j:(j + p$categories[p$categories[j, 4], 2] - 1))] <- matrix((model.matrix(~cat.columns - 1)[,
                                                                                                                    -1]), ncol = p$categories[p$categories[j, 4], 2], nrow = nrow(p$data))
              remove("cat.columns")
            }

            ## optional post-processing
            cmd <- p$post[j]  # SvB Aug 2009
            if (cmd != "") {
              eval(parse(text = cmd))
              p$data[!r[, j], j] <- imp[[j]][, i]
            }
          }  # end j loop
        }  # end i loop
        k2 <- k - from + 1
        for (j in 1:length(visitSequence)) {
          jj <- visitSequence[j]
          if (!is.factor(data[, jj])) {
            chainVar[j, k2, ] <- apply(imp[[jj]], 2, var)
            chainMean[j, k2, ] <- colMeans(as.matrix(imp[[jj]]))  ##pm 04/02
          }
          if (is.factor(data[, jj])) {
            for (mm in 1:m) {
              nc <- as.integer(factor(imp[[jj]][, mm], levels = levels(data[, jj])))
              chainVar[j, k2, mm] <- var(nc)
              chainMean[j, k2, mm] <- mean(nc)
            }
          }
        }
      }  # end iteration loop
      if (printFlag)
        cat("\n")
    }
  return(list(iteration = maxit, imp = imp, chainMean = chainMean, chainVar = chainVar))
}
environment(sampler2)<-environment(mice)

#' Function for multiple imputation by chained equations allowing for "lmm" imputation method.
#'
#' Modified version of \code{mice} function from the \code{mice} package to allow for "lmm" imputation method.
#' @export
#' @keywords internal
#' @import mice
#' @importFrom stats as.formula model.extract pnorm pt qt quantile rbinom rexp rnorm runif terms uniroot

mice2<-function (data, m = 5, method = vector("character", length = ncol(data)),
                 predictorMatrix = (1 - diag(1, ncol(data))), visitSequence = (1:ncol(data))[apply(is.na(data),
                                                                                                   2, any)], form = vector("character", length = ncol(data)),
                 post = vector("character", length = ncol(data)), defaultMethod = c("pmm",
                                                                                    "logreg", "polyreg", "polr"), maxit = 5, diagnostics = TRUE,
                 printFlag = TRUE, seed = NA, imputationMethod = NULL, defaultImputationMethod = NULL,
                 data.init = NULL, ...)
{
  check.visitSequence <- function(setup) {
    nmis <- setup$nmis
    nvar <- setup$nvar
    visitSequence <- setup$visitSequence
    if (!is.numeric(visitSequence)) {
      code <- pmatch(visitSequence, c("roman", "arabic",
                                      "monotone", "revmonotone"))
      if (!is.na(code) && code == 1)
        visitSequence <- (1:nvar)[nmis > 0]
      if (!is.na(code) && code == 2)
        visitSequence <- rev((1:nvar)[nmis > 0])
      if (!is.na(code) && code == 3)
        visitSequence <- order(nmis)[(sum(nmis == 0) +
                                        1):length(nmis)]
      if (!is.na(code) && code == 4)
        visitSequence <- rev(order(nmis)[(sum(nmis ==
                                                0) + 1):length(nmis)])
      if (is.na(code))
        stop("Argument visitSequence not recognized.\n")
    }
    if (all(nmis[visitSequence] == 0))
      stop(paste("No missing values found."))
    flags <- nmis == 0 & is.element(1:nvar, visitSequence)
    if (any(flags))
      visitSequence <- visitSequence[!flags]
    visitSequence <- visitSequence[visitSequence <= nvar]
    visitSequence <- visitSequence[visitSequence >= 1]
    if (length(visitSequence) == 0)
      stop(paste("No missing values found."))
    setup$visitSequence <- visitSequence
    return(setup)
  }
  check.predictorMatrix <- function(setup) {
    pred <- setup$predictorMatrix
    varnames <- setup$varnames
    nmis <- setup$nmis
    nvar <- setup$nvar
    vis <- setup$visitSequence
    method <- setup$method
    post <- setup$post
    if (!is.matrix(pred))
      stop("Argument 'predictorMatrix' not a matrix.")
    if (nvar != nrow(pred) | nvar != ncol(pred))
      stop(paste("The predictorMatrix has", nrow(pred),
                 "rows and", ncol(pred), "columns. Both should be",
                 nvar, "."))
    dimnames(pred) <- list(varnames, varnames)
    diag(pred) <- 0
    isclassvar <- apply(pred == -2, 2, any)
    for (j in 1:nvar) {
      if (method[j] == "" & any(pred[, j] != 0) & nmis[j] >
            0) {
        out <- varnames[j]
        updateLog(out = out)
        pred[, j] <- 0
        vis <- vis[vis != j]
        post[j] <- ""
        if (isclassvar[j])
          stop("Removed an incomplete class variable.")
      }
      if (nmis[j] == 0 & any(pred[j, ] != 0))
        pred[j, ] <- 0
    }
    setup$predictorMatrix <- pred
    setup$visitSequence <- vis
    setup$post <- post
    return(setup)
  }
  check.method <- function(setup, data) {
    method <- setup$method
    defaultMethod <- setup$defaultMethod
    visitSequence <- setup$visitSequence
    nmis <- setup$nmis
    nvar <- setup$nvar
    if (all(method == "")) {
      for (j in visitSequence) {
        y <- data[, j]
        if (is.numeric(y))
          method[j] <- defaultMethod[1]
        else if (nlevels(y) == 2)
          method[j] <- defaultMethod[2]
        else if (is.ordered(y) & nlevels(y) > 2)
          method[j] <- defaultMethod[4]
        else if (nlevels(y) > 2)
          method[j] <- defaultMethod[3]
        else if (is.logical(y))
          method[j] <- defaultMethod[2]
        else method[j] <- defaultMethod[1]
      }
    }
    if (length(method) == 1) {
      if (is.passive(method))
        stop("Cannot have a passive imputation method for every column.")
      method <- rep(method, nvar)
    }
    if (length(method) != nvar)
      stop(paste("The length of method (", length(method),
                 ") does not match the number of columns in the data (",
                 nvar, ").", sep = ""))
    active <- !is.passive(method) & nmis > 0 & !(method ==
                                                   "")
    passive.check <- is.passive(method) & nmis > 0 & !(method ==
                                                         "")
    check <- all(active == FALSE) & any(passive.check !=
                                          FALSE)
    if (check)
      fullNames <- rep("mice.impute.passive", length(method[passive.check]))
    else fullNames <- paste("mice.impute", method[active],
                            sep = ".")
    notFound <- !sapply(fullNames, exists, mode = "function",
                        inherit = TRUE)
    if (any(notFound))
      stop(paste("The following functions were not found:",
                 paste(fullNames[notFound], collapse = ", ")))
    for (j in visitSequence) {
      y <- data[, j]
      vname <- dimnames(data)[[2]][j]
      mj <- method[j]
      mlist <- list(m1 = c("logreg", "logreg.boot", "polyreg",
                           "lda", "polr"), m2 = c("norm", "norm.nob", "norm.predict",
                                                  "norm.boot", "mean", "2l.norm", "2L.norm", "2l.pan",
                                                  "2L.pan", "2lonly.pan", "quadratic", "ri"), m3 = c("norm",
                                                                                                     "norm.nob", "norm.predict", "norm.boot", "mean",
                                                                                                     "2l.norm", "2L.norm", "2l.pan", "2L.pan", "2lonly.pan",
                                                                                                     "quadratic", "logreg", "logreg.boot"))
      if (is.numeric(y) & (mj %in% mlist$m1))
        warning("Type mismatch for variable ", vname,
                "\nImputation method ", mj, " is for categorical data.",
                "\nIf you want that, turn variable ", vname,
                " into a factor,", "\nand store your data in a data frame.",
                call. = FALSE)
      else if (is.factor(y) & nlevels(y) == 2 & (mj %in%
                                                   mlist$m2))
        warning("Type mismatch for variable ", vname,
                "\nImputation method ", mj, " is not for factors.",
                call. = FALSE)
      else if (is.factor(y) & nlevels(y) > 2 & (mj %in%
                                                  mlist$m3))
        warning("Type mismatch for variable ", vname,
                "\nImputation method ", mj, " is not for factors with three or more levels.",
                call. = FALSE)
    }
    setup$method <- method
    return(setup)
  }
  check.data <- function(setup, data, allow.na = FALSE, ...) {
    pred <- setup$predictorMatrix
    nvar <- setup$nvar
    varnames <- setup$varnames
    meth <- setup$method
    vis <- setup$visitSequence
    isclassvar <- apply(pred == -2, 2, any)
    for (j in 1:nvar) {
      if (isclassvar[j] & is.factor(data[, j]))
        stop(paste("Class variable (column ", j, ") cannot be factor. Convert to numeric by as.integer()",
                   sep = ""))
    }
    for (j in 1:nvar) {
      if (!is.passive(meth[j])) {
            d.j <- data[, j]
            v <- if (is.character(d.j))
              NA
            else var(as.numeric(d.j), na.rm = TRUE)
            constant <- if (allow.na) {
              if (is.na(v))
                FALSE
              else (v < 1000 * .Machine$double.eps)
            }
            else is.na(v) || v < 1000 * .Machine$double.eps
        didlog <- FALSE
        if (constant & any(pred[, j] != 0)) {
          out <- varnames[j]
          pred[, j] <- 0
          updateLog(out = out, meth = "constant")
          didlog <- TRUE
        }
        if (constant & meth[j] != "") {
          out <- varnames[j]
          pred[j, ] <- 0
          if (!didlog)
            updateLog(out = out, meth = "constant")
          meth[j] <- ""
          vis <- vis[vis != j]
          post[j] <- ""
        }
      }
    }

    ispredictor <- apply(pred != 0, 2, any)
    if (any(ispredictor))
      droplist <- find.collinear(data[, ispredictor, drop = FALSE],
                                 ...)
    else droplist <- NULL
    if (length(droplist) > 0) {
      for (k in 1:length(droplist)) {
        j <- which(varnames %in% droplist[k])
        didlog <- FALSE
        if (any(pred[, j] != 0)) {
          out <- varnames[j]
          pred[, j] <- 0
          updateLog(out = out, meth = "collinear")
          didlog <- TRUE
        }
        if (meth[j] != "") {
          out <- varnames[j]
          pred[j, ] <- 0
          if (!didlog)
            updateLog(out = out, meth = "collinear")
          meth[j] <- ""
          vis <- vis[vis != j]
          post[j] <- ""
        }
      }
    }
    setup$predictorMatrix <- pred
    setup$visitSequence <- vis
    setup$post <- post
    setup$meth <- meth
    return(setup)
  }
  call <- match.call()
  if (!is.na(seed))
    set.seed(seed)
  if (!(is.matrix(data) | is.data.frame(data)))
    stop("Data should be a matrix or data frame")
  if ((nvar <- ncol(data)) < 2)
    stop("Data should contain at least two columns")
  data <- as.data.frame(data)
  nmis <- apply(is.na(data), 2, sum)
  if (sum(nmis) == 0)
    stop("No missing values found")
  varnames <- dimnames(data)[[2]]
  state <- list(it = 0, im = 0, co = 0, dep = "", meth = "",
                log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, co = 0, dep = "",
                             meth = "", out = "")
  if (!is.null(imputationMethod))
    method <- imputationMethod
  if (!is.null(defaultImputationMethod))
    defaultMethod <- defaultImputationMethod
  setup <- list(visitSequence = visitSequence, method = method,
                defaultMethod = defaultMethod, predictorMatrix = predictorMatrix,
                form = form, post = post, nvar = nvar, nmis = nmis, varnames = varnames)
  setup <- check.visitSequence(setup)
  setup <- check.method(setup, data)
  setup <- check.predictorMatrix(setup)
  setup <- check.data(setup, data, ...)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  p <- padModel(data, method, predictorMatrix, visitSequence,
                form, post, nmis, nvar)
  if (sum(duplicated(names(p$data))) > 0)
    stop("Column names of padded data should be unique")
  r <- (!is.na(p$data))
  imp <- vector("list", ncol(p$data))
  if (m > 0) {
    for (j in visitSequence) {
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(!r[,
                                                         j]), ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[r[, j] ==
                                                   FALSE], 1:m)
      y <- data[, j]
      ry <- r[, j]
      if (method[j] != "") {
        for (i in 1:m) {
          if (nmis[j] < nrow(data)) {
            if (is.null(data.init)) {
              imp[[j]][, i] <- mice.impute.sample(y,
                                                  ry, ...)
            }
            else {
              imp[[j]][, i] <- data.init[!ry, j]
            }
          }
          else imp[[j]][, i] <- rnorm(nrow(data))
        }
      }
    }
  }
  from <- 1
  to <- from + maxit - 1
  q <- survtd::sampler2(p, data, m, imp, r, visitSequence, c(from, to),
                printFlag, ...)
  for (j in p$visitSequence) p$data[(!r[, j]), j] <- NA
  imp <- q$imp[1:nvar]
  names(imp) <- varnames
  names(method) <- varnames
  names(form) <- varnames
  names(post) <- varnames
  names(visitSequence) <- varnames[visitSequence]
  if (!state$log)
    loggedEvents <- NULL
  if (state$log)
    row.names(loggedEvents) <- 1:nrow(loggedEvents)
  midsobj <- list(call = call, data = as.data.frame(p$data[,
                                                           1:nvar]), m = m, nmis = nmis, imp = imp, method = method,
                  predictorMatrix = predictorMatrix, visitSequence = visitSequence,
                  form = form, post = post, seed = seed, iteration = q$iteration,
                  lastSeedValue = .Random.seed, chainMean = q$chainMean,
                  chainVar = q$chainVar, loggedEvents = loggedEvents)
  if (diagnostics)
    midsobj <- c(midsobj, list(pad = p))
  oldClass(midsobj) <- "mids"
  return(midsobj)
}
environment(mice2)<-environment(mice)
