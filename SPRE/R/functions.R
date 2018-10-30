##' SPRE
##'
##' Perform SPRE analysis given predictors and responses
##' @title SPRE
##' @param predictor predictor vector in the linear model
##' @param response response vector in the linear model
##' @param info whether output fitting information
##' @return an object of class 'SPRE'
##' @author Deborah Weissman-Miller
##' @export
##' @importFrom stats lm
##' @importFrom stats pf
##' @references 
##' Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462.
##' Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184.
##' Weissman-Miller, D., & Graham, K. C. (2015). Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4(3), 161.
##' Weissman-Miller, D. (2016). On predicting survival in prostate cancer: using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE). International Journal of Statistics and Probability, 5(2), 19.
##' Weissman-Miller, D., Miller, R. J., & Shotwell, M. P. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational therapy international, 2017.
##' @examples
##' data('HEAT_stat')
##' modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
SPRE <- function(predictor, response, info = TRUE) {
    
    output_object <- list()
    class(output_object) <- "SPRE"
    
    output_object$predictor <- predictor
    output_object$response <- response
    
    lmfit <- lm(response ~ predictor)
    asummary <- summary(lmfit)
    sL <- asummary$coefficients[2, 1]
    sL <- abs(sL)
    
    output_object$sL <- sL
    output_object$lmfit <- lmfit
    
    fDist <- matrix(data = NA, nrow = length(response), ncol = 6, byrow = TRUE)
    colnames(fDist) <- c("Rsquare", "Fstatistic", "FstatNumer.", "FstatDenom.", "Pvalue", "Session")
    tempdata <- data.frame(predictor = predictor, response = response)
    
    for (i in length(response):1) {
        x <- summary(lm(response ~ predictor, data = tempdata))
        Rsq <- x$r.squared
        n <- x$fstatistic
        f <- as.numeric(x$fstatistic)
        P <- as.numeric(pf(f[1], f[2], f[3], lower.tail = F))
        SE <- predictor
        SE <- t(SE)
        tempdata <- tempdata[-c(nrow(tempdata)), ]
        fDist[i, ] <- c(Rsq, n, P, SE[i])
    }
    
    output_object$fDist <- fDist
    
    if (info) {
        cat("\n", "Step 2: Backwards stepwise OLS linear regression", "\n", "\n")
        print(signif(fDist, digits = 6))
        cat(" \n ")
    }
    
    tester <- apply(fDist, FUN = which.max, 2)
    y <- tester["Fstatistic"]
    q <- fDist[y, "Fstatistic"]
    qn <- fDist[y, "Session"]
    v <- fDist[y, "Session"]
    u <- fDist[y, "Session"] + 1
    R2 <- fDist[y, "Rsquare"]
    pv <- fDist[y, "Pvalue"]
    k <- abs(log(abs(1 - (sL * v))))
    
    output_object$y <- y
    output_object$k <- k
    output_object$q <- q
    output_object$qn <- qn
    output_object$pv <- pv
    output_object$R2 <- R2
    
    if (info) {
        cat("\n")
        cat("slope parameter: ")
        cat("betaHat1 =", sprintf("%.4f", sL), "\n")
        cat("shape parameter: ")
        cat("k_value -", sprintf("%.4f", k), "\n")
        cat("scale parameter is the value (tau) at the change point - given as the highest F statistic")
        cat("\n")
        cat("\n")
        cat("Session number for -", sprintf("%.4f", qn), "\n")
        cat("Session value for -", sprintf("%.4f", pv), "\n")
        cat("R^2 at the change point -", sprintf("%.4f", R2), "\n")
        cat("If R^2 is less than 0.45 at the change point, do not predict for inference")
        cat("\n")
    }
    
    r <- ((1 - exp(-(u/v)^k))/(1 - exp(-(v/v)^k))) == ((1 - exp(-(u/v)^k))/(1 - exp(-(v/v)^k)))
    
    Pred <- matrix(data = NA, nrow = 40, ncol = 5, byrow = TRUE)
    colnames(Pred) <- c("Session", "Ratio", "tau", "k.factor", "Predictions")
    
    stat <- response
    z <- response[y]
    
    output_object$z <- z
    
    for (i in 2:40) {
        if (stat[2] < stat[length(stat)]) {
            v <- i + fDist[y, "Session"]
            u <- i + (fDist[y, "Session"] + 1)
            sess <- i + fDist[y, "Session"]
            z[i] <- response[y] * r
        } else {
            if (stat[2] > stat[length(stat)]) 
                v <- i + (fDist[y, 6] + 1)
            u <- i + fDist[y, 6]
            sess <- i + fDist[y, 6]
            z[i] <- response[y] * r
        }
        r <- ((1 - exp(-(u/v)^k))/(1 - exp(-(v/v)^k)))
        z[i] <- z[i - 1] * r
        
        Pred[i, ] = c(sess, r, fDist[y, 6], k, z[i])
    }
    if (info) {
        cat("\n", "\n")
        print(signif(Pred, digits = 6))
        cat("\n")
        cat("Effective clinical stability may be determined when the second significant digit (after the decimal place) remains constant for 2 predictions.", 
            "These 2 sessions, or more, are then effective clinical stability as long as treatment is maintained to these sessions.")
        cat("\n", "\n")
    }
    
    output_object$Pred <- Pred
    
    return(output_object)
}

##' print SPRE
##'
##' print an SPRE object
##' @title print SPRE
##' @param x an object of class 'SPRE'
##' @param ... arguments passed to print.default
##' @return output some information to console
##' @author Deborah Weissman-Miller
##' @export
##' @references 
##' Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462.
##' Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184.
##' Weissman-Miller, D., & Graham, K. C. (2015). Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4(3), 161.
##' Weissman-Miller, D. (2016). On predicting survival in prostate cancer: using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE). International Journal of Statistics and Probability, 5(2), 19.
##' Weissman-Miller, D., Miller, R. J., & Shotwell, M. P. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational therapy international, 2017.

##' @examples
##' data('HEAT_stat')
##' modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
##' print(modelfit)
print.SPRE <- function(x, ...) {
    print(signif(x$fDist, digits = 6))
    print(signif(x$Pred, digits = 6))
}

##' plot residuals
##'
##' residuals plot from the linear model given in the first SPRE step
##' @title plot residuals
##' @param SPRE_object an object of class 'SPRE'
##' @return output some plots to device
##' @author Deborah Weissman-Miller
##' @export
##' @importFrom graphics plot
##' @importFrom graphics abline
##' @importFrom graphics par
##' @references 
##' Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462.
##' Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184.
##' Weissman-Miller, D., & Graham, K. C. (2015). Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4(3), 161.
##' Weissman-Miller, D. (2016). On predicting survival in prostate cancer: using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE). International Journal of Statistics and Probability, 5(2), 19.
##' Weissman-Miller, D., Miller, R. J., & Shotwell, M. P. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational therapy international, 2017.
##' @examples
##' data('HEAT_stat')
##' modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
##' plot_residuals(modelfit)
plot_residuals <- function(SPRE_object) {
    opar <- par("mfrow", "mar")
    par(mfrow = c(2, 1), mar = c(4.7, 5.8, 2.5, 3.8))
    
    plot(SPRE_object$lmfit, 1, ann = FALSE, pch = 18, title("SPRE Residuals", col.main = "blue", xlab = "Fitted Values", 
        ylab = "Residuals"))
    abline(0, 0)
    
    plot(SPRE_object$lmfit, 2, ann = FALSE, pch = 18, title("Normal Q-Q", col.main = "blue", xlab = "xlab
             Quantiles", 
        ylab = "Standardized Residuals"))
    abline(0, 0)
    
    par(opar)
}

##' plot predictions
##'
##' plot the predictions from SPRE model
##' @title plot predictions
##' @param SPRE_object an object of class 'SPRE'
##' @return output plots to device
##' @author Deborah Weissman-Miller
##' @export
##' @importFrom graphics plot
##' @importFrom graphics lines
##' @importFrom graphics legend
##' @importFrom graphics title
##' @references 
##' Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462.
##' Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184.
##' Weissman-Miller, D., & Graham, K. C. (2015). Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4(3), 161.
##' Weissman-Miller, D. (2016). On predicting survival in prostate cancer: using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE). International Journal of Statistics and Probability, 5(2), 19.
##' Weissman-Miller, D., Miller, R. J., & Shotwell, M. P. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational therapy international, 2017.
##' @examples
##' data('HEAT_stat')
##' modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
##' plot_predictions(modelfit)
plot_predictions <- function(SPRE_object) {
    plot(SPRE_object$Pred[, "Predictions"] ~ SPRE_object$Pred[, "Session"], ann = FALSE, type = "o", 
        pch = 19, col = "blue")
    lines(SPRE_object$Pred[, "Predictions"] ~ SPRE_object$Pred[, "Session"], lwd = 2, col = "blue")
    legend("bottomleft", legend = "ChangePt       ", cex = 0.8)
    title("SPRE Predictions from a 'change point'", col.main = "blue", xlab = "Prediction Sessions", 
        ylab = "Predictions")
}

##' stability analysis
##'
##' stability analysis of SPRE predictions
##' @title stability analysis
##' @param SPRE_object an object of class 'SPRE'
##' @return output information to console and plots to device
##' @author Deborah Weissman-Miller
##' @export
##' @importFrom graphics curve
##' @importFrom graphics par
##' @importFrom stats pweibull
##' @references 
##' Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462.
##' Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184.
##' Weissman-Miller, D., & Graham, K. C. (2015). Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4(3), 161.
##' Weissman-Miller, D. (2016). On predicting survival in prostate cancer: using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE). International Journal of Statistics and Probability, 5(2), 19.
##' Weissman-Miller, D., Miller, R. J., & Shotwell, M. P. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational therapy international, 2017.
##' @examples
##' data('HEAT_stat')
##' modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
##' stability(modelfit)
stability <- function(SPRE_object) {
    Pred2 <- matrix(SPRE_object$Pred[, "Predictions"])
    Pr <- format(round(Pred2, 2), nsmall = 2)
    cat("Predictions, significant digits = 2, rows 1:20:", "\n", "\n")
    print(Pr[1:20])
    
    dup <- format(round(Pred2, 2), nsmall = 2)
    xu1 <- dup[duplicated(dup, fromLast = TRUE)]
    xu2 <- xu1[1]
    cat("\n")
    cat("First predicted duplicate value:", xu2[1])
    
    dup2 <- which(dup[, 1] == xu2[1])
    cat("\n")
    cat("Row numbers for predicted sessions:", dup2)
    cat("\n")
    
    x <- SPRE_object$Pred[, "Predictions"]
    
    opar <- par("mfrow", "mar")
    par(mfrow = c(1, 1), mar = c(5.5, 5.2, 5.5, 3))
    m2 <- SPRE_object$qn
    pv <- SPRE_object$pv
    k <- SPRE_object$k
    sL <- SPRE_object$sL
    plot(x, main = "Probability for SPRE Predictions\nwith\neffective clinical stability (ECS)", xlab = paste("Prediction & Pvalue from Change Point =", 
        m2, ", p =", sprintf("%.4f", pv)), ylab = "Probability", xlim = c(0, 40), ylim = c(0, 1))
    
    cat("\n")
    cat("Clinical stability is determined from the probability curve for SPRE predictions,", "if any 2 predictions from the second significant digit (after the decimal place) remains constant.", 
        "Otherwise, there is clinical significance at the vertical line (zero) at the change point.")
    
    curve(pweibull(x, scale = sL, shape = k), from = 2, to = 40, add = TRUE, lwd = 2)
    
    cat("\n", "\n")
    cat("First six predicted probability values", "\n")
    
    m4 <- pweibull(x, shape = k, scale = sL, lower.tail = TRUE, log.p = FALSE)[1:6]
    print(pweibull(x, shape = k, scale = sL, lower.tail = TRUE, log.p = FALSE)[1:6])
    abline(v = 0, col = "blue", lty = 3)
    m3 <- dup2
    abline(v = m3, col = "red", lty = 2, lwd = 1.5)
    par(opar)
}

##' plot weibull
##'
##' plot the weibull distribution from SPRE model
##' @title plot weibull
##' @param SPRE_object an object of class 'SPRE'
##' @return output plots to device
##' @author Deborah Weissman-Miller
##' @export
##' @references 
##' Weissman-Miller, D., Shotwell, M.P. & Miller, R.J. (2012).  New single-subject and small-n design in occupational therapy:  Application to weight loss in obesity.  American Journal of Occupational Therapy, 66, 455-462.
##' Weissman-Miller, D. (2013). Novel point estimation from a Semiparametric Ratio Estimator (SPRE): Long-term health outcomes from short-term linear data, with application to weight loss in obesity. International Journal of Biostatistics, 9(2): 175-184.
##' Weissman-Miller, D., & Graham, K. C. (2015). Novel scale development for fear of falling and falls: Analysed using a Semiparametric Ratio Estimator (SPRE). International Journal of Statistics and Probability, 4(3), 161.
##' Weissman-Miller, D. (2016). On predicting survival in prostate cancer: using an extended maximum spacing method at the change point of the semiparametric ratio estimator (SPRE). International Journal of Statistics and Probability, 5(2), 19.
##' Weissman-Miller, D., Miller, R. J., & Shotwell, M. P. (2017). Translational Research for Occupational Therapy: Using SPRE in Hippotherapy for Children with Developmental Disabilities. Occupational therapy international, 2017.
##' @examples
##' data('HEAT_stat')
##' modelfit<-SPRE(predictor=HEAT_stat$Session,response=HEAT_stat$FData)
##' plot_weibull(modelfit)
plot_weibull <- function(SPRE_object) {
    z <- sort(SPRE_object$Pred[-1, "Predictions"])
    fit.weibull <- fitdistrplus::fitdist(z, "weibull", lower = c(0, 0))
    plot(fit.weibull, cex.main = 0.8)
}
