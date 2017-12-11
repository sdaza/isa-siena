# extract function for metafor

extract.rma <- function(model, include.nobs = TRUE, aic = TRUE, bic = TRUE, tau = TRUE,
	cochran = TRUE, i2 = FALSE, ...) {
  
  names <- rownames(model$b)
  co <- model$b[, 1]
  se <- model$se
  pval <- model$pval
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (aic == TRUE) {
    maic <- AIC.rma(model)
    gof <- c(gof, maic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  if (bic == TRUE) {
    mbic <- BIC.rma(model)
    gof <- c(gof, mbic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, FALSE)
  }
    if (tau == TRUE) {
    i2 <- model$tau2
    gof <- c(gof, i2)
    gof.names <- c(gof.names, "$\\tau^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (cochran == TRUE) {
    QEp <- model$QEp 
    gof <- c(gof, QEp)
    gof.names <- c(gof.names, "Cochran Q-test (p-value)")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (i2 == TRUE) {
    i2 <- model$I2
    gof <- c(gof, i2)
    gof.names <- c(gof.names, "$I^2$ statistic")
    gof.decimal <- c(gof.decimal, FALSE)
  }
if ( include.nobs == TRUE) {
  n <- model$k
  gof <- c(gof, n)
  gof.names <- c(gof.names, "Num.\\ obs.")
  gof.decimal <- c(gof.decimal, FALSE)
}
  
  tr <- createTexreg(  
    coef.names = names, 
    coef = co, 
    se = se, 
    pvalues = pval, 
    gof.names = gof.names, 
    gof = gof, 
    gof.decimal = gof.decimal
    )

  return(tr)
}

setMethod("extract", signature = className("rma.uni", "rma"), 
          definition = extract.rma)