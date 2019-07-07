anova.ltm <-
    function (object, object2, conditionalLR = FALSE, ...) {
        if (!inherits(object, "ltm"))
            stop("Use only with 'ltm' objects.\n")
        if (missing(object2))
            stop("anova.ltm() computes LRTs between two nested IRT models.")
        if (!inherits(object2, "ltm") && !inherits(object2, "tpm"))
            stop("'object2' must inherit from either class 'ltm' or class 'tpm'.\n")
        if (!isTRUE(all.equal(object$X, object2$X)))
            warning("it seems that the two models are fitted in different data sets.")
        L0 <- logLik(object)
        L1 <- logLik(object2)
        nb0 <- attr(L0, "df")
        nb1 <- attr(L1, "df")
        df. <- nb1 - nb0
        if (df. < 0)
            stop("'object' is not nested in 'object2'.\n")
        LRT <- -2 * (L0 - L1)
        attributes(LRT) <- NULL
        if (LRT < 0)
            warning( "either the two models are not nested or the model represented by 'object2' fell on a local maxima.\n" )
        p.value <- pchisq(LRT, df., lower.tail = FALSE)
        
        out <-
            list(
                nam0 = deparse(substitute(object)),
                L0 = L0,
                aic0 = AIC(object),
                bic0 = AIC(object, k = log(attr(L0, "n"))),
                nam1 = deparse(substitute(object2)),
                L1 = L1,
                aic1 = AIC(object2),
                bic1 = AIC(object2, k = log(attr(L1, "n"))),
                LRT = LRT,
                df = df.,
                p.value = p.value,
                call = object$call
            )
        class(out) <- "aov.ltm"
        
        if (conditionalLR == FALSE) {
            return(out)
        }
        if (conditionalLR == TRUE) {
            cs <- coef(object2)[, 1]
            sens <- c(as.numeric(as.character(paste( "1e-", 1:8, sep = "" ))), 0)
            p_values <- data.frame("Sensibility" = sens)
            if(length(object2$constraint != 0)){
                warning("When using constrained models we advise to check that the degrees of freedom have been calculated correctly.\n")
            }
            constraints <- length(object2$constraint[,1])
            for (i in 1:length(sens)) {
                dof <- df. - sum(cs < sens[i]) + length(object2$constraint[object2$constraint[,2]==1 ,1])
                if (dof < 0){
                    stop( "Something wrong has happened.\n This implementation is still under development, please fill an issue on github (cobac/ltm) with the steps required to reproduce this error.\n" )
                    }
                p_values$p.values[i] <- round(pchisq(LRT, dof, lower.tail = FALSE), 3)
                p_values$df[i] <- dof
            }
            return(
                list(
                    "Naive Likelihood Ratio Test" = out,
                    "Conditional Likelihood Ratio Test" = p_values
                )
            )
        }
    }
