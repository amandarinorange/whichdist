####################################################
#### whichdist
####################################################
#' @name whichdist
#' @aliases whichdist
#' @title Assess the Best-Fitting Distribution For a Continuous Variable
#' @description With this function, the user can input a continuous count metric they aim to use as the outcome variable in an analysis and assess which type of regression to run
#' @usage whichdist(df, countvar, diagnostics = FALSE)
#' @param df The dataframe in which the count variable is a column
#' @param countvar - A column of integer values -- i.e., a count metric -- that the user intends to evaluate
#' @param diagnostics - Boolean, whether the function should return goodness-of-fit metrics in addition to the recommended regression model. Default is FALSE.
#' @return A list object with the recommended regression and best-fitting statistical distribution. Also returns results of goodness of fit and estimated parameters when diagnostics is set to TRUE.
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom MASS "glm.nb"
#' @importFrom pscl "zeroinfl"
#' @importFrom fitdistrplus "fitdist"
#' @importFrom stats "pchisq"
#' @importFrom pscl "zeroinfl"
#' @importFrom gamlss.dist "ZIP"
#' @importFrom emdbook "dzinbinom"


whichdist <- function(df, countvar, diagnostics = FALSE) {
  suppressWarnings({
    options(stringsAsFactors = F)

    col_ind <- which(colnames(df) == countvar)

    # evaluate whether the metric is an integer; if yes, then apply count-specific distributions
    if (sum(df[!is.na(df[,col_ind]),col_ind] - floor(df[!is.na(df[,col_ind]),col_ind])!=0)!=0) {
      return("This function only works for count (integer) variables")
    }
    else if (sum(is.na(df[,col_ind])) > 0) {
      return("Count variable cannot contain missing values")
    }
    else if (sum(df[!is.na(df[,col_ind]),col_ind] - floor(df[!is.na(df[,col_ind]),col_ind])!=0)==0) {
      glmnbfit <- glm.nb(df[,col_ind] ~ 1)
      # make output table
      headers.1 <- c("normal_fitdistr","pois_fitdistr", "nbinom_fitdistr", "zifpois_fitdistr", "zifnbinom_fitdistr")
      metric <- c("loglikelihood", "aic")

      gof.res <- cbind(metric, data.frame(matrix(vector(), nrow = length(metric),length(headers.1),
                                                 dimnames=list(c(metric),c(headers.1)))))

      headers_lrt <- c('LRT_PvsNorm','LRT_NBvsNorm',
                       'LRT_NBvsP', 'LRT_ZINBvsNB',
                       'LRT_ZIPvsP', 'LRT_ZINBvsZIP')
      lrt_res <- data.frame(matrix(vector(), nrow = 1, ncol = length(headers_lrt), dimnames = list("LRT p-value", headers_lrt)))

      headersests <- c("normal_mean","normal_sd",
                       "pois_lambda", "pois_sd",
                       "nbinom_size", "nbinom_mu",
                       "zifpois_mu", "zifpois_sigma",
                       "zifnbinom_mu", "zifnbinom_sigma", "zifnbinom_nu")
      estimate.res <- data.frame(matrix(vector(), 1, length(headersests),dimnames=list(c(countvar),c(headersests))))

      # zero measurement indicator
      counts0 <- df[,col_ind] == 0
      # number of samples
      nsamples=nrow(df)
      # number of nonzero samples
      nn0 = sum(df[,col_ind]>0)
      #dropout
      p0 = mean(df[,col_ind]==0)
      # mean and dispersion of nonzero portion of normalised counts
      estmu = sum((!counts0) * df[,col_ind])/nn0
      s2 = sum((!counts0) * (df[,col_ind] - estmu)^2)/(nn0 - 1)
      disp =  1 / (estmu^2/(s2 - estmu + 1e-04))

      # fit lm regression with fitdistrplus
      tmp.fit.lm.wo <- try(fitdistrplus::fitdist(data = df[,col_ind], 'norm'), silent = TRUE)

      # fit poisson with fitdistrplus
      tmp.fit.pois.wo <- try(fitdistrplus::fitdist(data = df[,col_ind], 'pois'), silent = TRUE)

      # fit negative binomial with fitdistrplus
      tmp.fit.nbinom.wo <- try(fitdistrplus::fitdist(data = df[,col_ind], 'nbinom'), silent = TRUE)

      if (min(df[,col_ind])==0) {
        # fit zero inflated poisson with fitdistrplus
        tmp.fit.zifpois.wo <- try(fitdistrplus::fitdist(df[,col_ind], "ZIP", start = list(mu = estmu, sigma=p0),
                                                        discrete = TRUE, lower = c(0, 0), upper = c(Inf, 1)), silent = TRUE)
        # fit zero inflated negative binomial with fitdistrplus
        tmp.fit.zifnbinom.wo <- try(fitdistrplus::fitdist(df[,col_ind], "ZINBI", start = list(mu = estmu, sigma=disp, nu=p0),
                                                          discrete = TRUE, lower = c(0, 0, 0), upper = c(Inf,Inf, 1)), silent = TRUE)
      }


      # lm with fitdistr
      lm.gof.fitdistr <- summary(tmp.fit.lm.wo)
      lmfitdistr.loglik <-  try(lm.gof.fitdistr$loglik)
      lmfitdistr.aic <-  try(lm.gof.fitdistr$aic)
      lmfitdistr_est <- try(tmp.fit.lm.wo$estimate)
      gof.res[, grep("^normal_fitdistr", colnames(gof.res), perl=TRUE)] <-  rbind(lmfitdistr.loglik, lmfitdistr.aic)
      estimate.res[, grep("^normal", colnames(estimate.res), perl=TRUE)] <- rbind(lmfitdistr_est[1],lmfitdistr_est[2])

      # poisson with fitdistr
      pois.gof.fitdistr <- summary(tmp.fit.pois.wo)
      poisfitdistr.loglik <-  try(pois.gof.fitdistr$loglik)
      poisfitdistr.aic <- try(pois.gof.fitdistr$aic)
      poisfitdistr_est <- try(tmp.fit.pois.wo$estimate)
      poisfitdistr_sd <- try(tmp.fit.pois.wo$sd)
      gof.res[, grep("^pois_fitdistr", colnames(gof.res), perl=TRUE)] <-  rbind(poisfitdistr.loglik, poisfitdistr.aic)
      estimate.res[, grep("^pois", colnames(estimate.res), perl=TRUE)] <- rbind(poisfitdistr_est, poisfitdistr_sd)

      # neg binom with fitdistr
      nbinom.gof.fitdistr <- try(summary(tmp.fit.nbinom.wo))
      nbinomfitdistr.loglik <-  try(nbinom.gof.fitdistr$loglik)
      nbinomfitdistr.aic <- try(nbinom.gof.fitdistr$aic)
      nbinomfitdistr_est <- try(tmp.fit.nbinom.wo$estimate)
      gof.res[, grep("^nbinom_fitdistr", colnames(gof.res), perl=TRUE)] <-  rbind(nbinomfitdistr.loglik, nbinomfitdistr.aic)
      estimate.res[, grep("^nbinom", colnames(estimate.res), perl=TRUE)] <- nbinomfitdistr_est


      if (min(df[, col_ind])==0) {
        # zero inflated poisson with fitdistr
        zifpois.gof.fitdistr <- summary(tmp.fit.zifpois.wo)
        zifpoisfitdistr.loglik <-  zifpois.gof.fitdistr$loglik
        zifpoisfitdistr.aic <-  zifpois.gof.fitdistr$aic
        zifpoisfitdistr_est <- tmp.fit.zifpois.wo$estimate
        gof.res[, grep("^zifpois_fitdistr", colnames(gof.res), perl = TRUE)] <- rbind(zifpoisfitdistr.loglik, zifpoisfitdistr.aic)
        estimate.res[, grep("^zifpois", colnames(estimate.res), perl = TRUE)] <- zifpoisfitdistr_est

        # zero inflated neg binom with fitdistr
        zifnbinom.gof.fitdistr <- summary(tmp.fit.zifnbinom.wo)
        zifnbfitdistr.loglik <-  zifnbinom.gof.fitdistr$loglik
        zifnbfitdistr.aic <-  zifnbinom.gof.fitdistr$aic
        zifnbfitdistr_est <- tmp.fit.zifnbinom.wo$estimate
        gof.res[, grep("^zifnbinom_fitdistr", colnames(gof.res), perl = TRUE)] <- rbind(zifnbfitdistr.loglik, zifnbfitdistr.aic)
        estimate.res[, grep("^zifnb", colnames(estimate.res), perl = TRUE)] <- zifnbfitdistr_est
      }
      # LR Test for P > Norm ?
      tmp.lrt.pnorm = try(pchisq(as.numeric(2 * (logLik(tmp.fit.pois.wo ) - logLik(tmp.fit.lm.wo))),
                                 df = length(tmp.fit.pois.wo$estimate)-length(tmp.fit.lm.wo$estimate),
                                 lower.tail = FALSE))
      lrt_res[, grep("LRT_PvsNorm", colnames(lrt_res), perl=TRUE)] <- tmp.lrt.pnorm

      # LR Test for NB > Norm ?
      tmp.lrt.nbnorm = try(pchisq(as.numeric(2 * (logLik(tmp.fit.nbinom.wo ) - logLik(tmp.fit.lm.wo))),
                                  df = length(tmp.fit.nbinom.wo$estimate)-length(tmp.fit.lm.wo$estimate),
                                  lower.tail = FALSE))
      lrt_res[, grep("LRT_NBvsNorm", colnames(lrt_res), perl=TRUE)] <- tmp.lrt.nbnorm

      # LR Test for NB > P ?
      tmp.lrt.nbp = try(pchisq(as.numeric(2 * (logLik(tmp.fit.nbinom.wo) - logLik(tmp.fit.pois.wo))),
                               df = length(tmp.fit.nbinom.wo$estimate)-length(tmp.fit.pois.wo$estimate),
                               lower.tail = FALSE))
      lrt_res[, grep("LRT_NBvsP", colnames(lrt_res), perl=TRUE)] <- tmp.lrt.nbp

      if (min(df[, col_ind])==0) {
        # LR Test for ZINB > NB ?
        tmp.lrt.nbzinb = try(pchisq(as.numeric(2 * (logLik(tmp.fit.zifnbinom.wo) - logLik(tmp.fit.nbinom.wo))),
                                    df = length(tmp.fit.zifnbinom.wo$estimate)-length(tmp.fit.nbinom.wo$estimate), lower.tail = FALSE))
        lrt_res[, grep("LRT_ZINBvsNB", colnames(lrt_res), perl=TRUE)] <- tmp.lrt.nbzinb

        # LR Test for ZIP > P ?
        tmp.lrt.zipp = try(pchisq(as.numeric(2 * (logLik(tmp.fit.zifpois.wo) - logLik(tmp.fit.pois.wo))),
                                  df = length(tmp.fit.zifpois.wo$estimate)-length(tmp.fit.pois.wo$estimate),
                                  lower.tail = FALSE))
        lrt_res[, grep("LRT_ZIPvsP", colnames(lrt_res), perl=TRUE)] <- tmp.lrt.zipp

        # LR Test for ZINB > ZIP ?
        tmp.lrt.zinbzip = try(pchisq(as.numeric(2 * (logLik(tmp.fit.zifnbinom.wo) - logLik(tmp.fit.zifpois.wo))),
                                     df = length(tmp.fit.zifnbinom.wo$estimate)-length(tmp.fit.zifpois.wo$estimate),
                                     lower.tail = FALSE))
        lrt_res[, grep("LRT_ZINBvsZIP", colnames(lrt_res), perl=TRUE)] <- tmp.lrt.zinbzip
      }

      lowest_aic_model <- colnames(gof.res[gof.res$metric == "aic",c(2:ncol(gof.res))])[apply(gof.res[gof.res$metric == "aic",c(2:ncol(gof.res))], 1, which.min)]
      if (grepl(pattern = "normal", x = lowest_aic_model, ignore.case = T)) {
        suggested_model <- "normal distribution: try a t-test, linear regression, or ANOVA"
      }

      if (grepl(pattern = "^pois", x = lowest_aic_model, ignore.case = T)) {
        suggested_model <- "Poisson distribution: try glm() with family = 'poisson'"
      }

      if (grepl(pattern = "^nbinom", x = lowest_aic_model, ignore.case = T)) {
        suggested_model <- "negative binomial distribution: try glm.nb() from MASS package"
      }

      if (grepl(pattern = "^zifpois", x = lowest_aic_model, ignore.case = T)) {
        suggested_model <- "zero-inflated Poisson distribution: try zeroinfl() from pscl package with dist = 'poisson'"
      }

      if (grepl(pattern = "^zifnbinom", x = lowest_aic_model, ignore.case = T)) {
        suggested_model <- "zero-inflated negative binomial distribution: try zeroinfl() from pscl package with dist = 'negbin'"
      }

      # list result object:
      if (diagnostics == TRUE) {
        return(list(goodness_of_fit_metrics=gof.res, likelihood_ratio_test_results = lrt_res, estimates=estimate.res,
                    paste0("results suggest that your data best fit a ", suggested_model)))
      }
      if (diagnostics == FALSE) {
        return(paste0("results suggest that your data best fit a ", suggested_model))
      }
    }
  })
}

