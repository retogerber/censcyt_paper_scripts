#' confidence interval from a mice::mipo object
#'
#' @param pooled_fit \code{\link[mice]{mipo}} object
#' @param alpha cutoff
#'
#' @return data.frame with lower and upper ci values
#' @export
#'
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, lm_formula, type = "lm", n_levels_fixeff=2)
#'  cmi_out <- conditional_multiple_imputation(data,lm_formula)
#'  mipo <- mice::pool(cmi_out$fits)
#'  ci_from_mipo(mipo)
ci_from_mipo <- function(pooled_fit, alpha=0.05){
  stopifnot(mice::is.mipo(pooled_fit))
  df <- pooled_fit$pooled$df
  estimate <- pooled_fit$pooled$estimate
  t <- pooled_fit$pooled$t
  data.frame(ci_lower = estimate - qt(p = 1-alpha/2,df = df)*sqrt(t),
             ci_upper = estimate + qt(p = 1-alpha/2,df = df)*sqrt(t),
             row.names = rownames(pooled_fit$pooled))
}


#' p-values and ci from fits
#'
#' p-values and confidence intervals from a list of fits. Additionally various
#' measures from combining of multiple imputation are calculated.
#'
#' @param fit_list list of fits
#'
#' @return data.frame
#' @export
#'
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, lm_formula, type = "lm", n_levels_fixeff=2)
#'  cmi_out <- conditional_multiple_imputation(data,lm_formula)
#'  pval_and_ci_from_fits(cmi_out$fits)
pval_and_ci_from_fits <- function(fit_list){
  if (!is.null(fit_list) & is.list(fit_list) & (length(fit_list)>1) & !(class(fit_list)[1] %in% c("lm","glm"))){
    tmp_mipo <- mice::pool(fit_list)
    pval <- summary(tmp_mipo)[ ,"p.value"]
    names(pval) <-  rownames(coef(summary(fit_list[[1]])))
    ci <- ci_from_mipo(tmp_mipo)
    estimate <- tmp_mipo$pooled[ ,"estimate"]
    total_var <- tmp_mipo$pooled[ ,"t"]
    within_var <- tmp_mipo$pooled[ ,"ubar"]
    between_var <- tmp_mipo$pooled[ ,"b"]
    df_adj <- tmp_mipo$pooled[ ,"df"]
    riv <- tmp_mipo$pooled[ ,"riv"]
    lambda <- tmp_mipo$pooled[ ,"lambda"]
    fmi <- tmp_mipo$pooled[ ,"fmi"]
    std_err <- summary(tmp_mipo)[ ,"std.error"]
  } else if (length(fit_list) == 1  | (class(fit_list)[1] %in% c("lm","glm"))) {
    pval <- coef(summary(fit_list))[ ,4]
    ci <- suppressMessages(confint(fit_list, method = "Wald"))
    estimate <-  coef(summary(fit_list))[ ,1]
    std_err <- coef(summary(fit_list))[ ,2]
  }
  return(
    data.frame(
      estimate = estimate,
      std_err =  tryCatch(std_err, error = function(e) rep(NA,length(pval))),
      total_var = tryCatch(total_var, error = function(e) rep(NA,length(pval))),
      within_var = tryCatch(within_var, error = function(e) rep(NA,length(pval))),
      between_var = tryCatch(between_var, error = function(e) rep(NA,length(pval))),
      df_adj = tryCatch(df_adj, error = function(e) rep(NA,length(pval))),
      riv = tryCatch(riv, error = function(e) rep(NA,length(pval))),
      lambda = tryCatch(lambda, error = function(e) rep(NA,length(pval))),
      fmi = tryCatch(fmi, error = function(e) rep(NA,length(pval))),
      ci_lower = tryCatch(ci[seq_along(pval), 1], error = function(e) rep(NA,length(pval))),
      ci_upper = tryCatch(ci[seq_along(pval), 2], error = function(e) rep(NA,length(pval))),
      pval = tryCatch(pval, error = function(e) rep(NA,length(pval))),
      row.names = tryCatch(names(pval), error = function(e) rep(NA,length(pval)))
    )
  )
}

#' Log ratio of censoring
#'
#' Log ratio of censoring dependent on a covariate with two levels
#'
#' @param data data.frame
#' @param covariate name of the censored covariate, rowname of 'data'
#' @param censoring_indicator name of censoring indicator, rowname of 'data'
#'
#' @return double
#' @export
#'
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, lm_formula, type = "lm",n_levels_fixeff=2)
#'  log_ratio_censoring_per_covariate_level(data,"Z","I")
log_ratio_censoring_per_covariate_level <- function(data, covariate, censoring_indicator){
  unique_covs <- unique(data[[covariate]])
  lens <- purrr::map(1:2, ~ dim(data[data[[covariate]] == unique_covs[.x], ])[1])
  n_cens <- purrr::map(1:2, ~ sum(data[[censoring_indicator]][data[[covariate]] == unique_covs[.x]]))
  numerator <- max((lens[[1]]-n_cens[[1]]),1e-10)
  denominator <- max((lens[[2]]-n_cens[[2]]),1e-10)
  return(log(numerator/denominator))
}

#' MSE of fits
#'
#' mean squared error of the response for a list of fits
#'
#' @param fit_list list of fits
#' @param response name of the response
#'
#' @return double
#' @export
#'
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, lm_formula, type = "lm", n_levels_fixeff=2)
#'  cmi_out <- conditional_multiple_imputation(data,lm_formula)
#'  mean_mse_of_fits(cmi_out$fits,"Y")
mean_mse_of_fits <- function(fit_list,response){
  if (!is.null(fit_list) & is.list(fit_list) & (length(fit_list)>1) & !(class(fit_list)[1] %in% c("lm","glm"))){
    purrr::map(fit_list, function(fit) {
      response_values <- model.frame(fit)[[response]]
      sum((response_values - fitted.values(fit)) ^ 2) / length(response_values)
    }) %>%
      unlist() %>%
      mean() %>%
      return()
  } else if (length(fit_list) == 1  | (class(fit_list)[1] %in% c("lm","glm"))) {
    response_values <- model.frame(fit_list)[[response]]
    return(sum((response_values - fitted.values(fit_list)) ^ 2) / length(response_values))
  } else {
    return(NA)
  }
}

#' MSE
#'
#' Mean squared error of response
#'
#' @param data data.frame
#' @param response column name to calculate mse from
#'
#' @return double
#' @export
#'
#' @examples
#'  lm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, lm_formula, type = "lm")
#'  mse_of_data_sim(data,"Y")
mse_of_data_sim <- function(data, response){
  sum((data[[response]] - data[[paste0(response,"_True")]])^2)/dim(data)[1]
}

#' General signal to noise ratio
#'
#' General signal to noise ratio for generalized linear models
#'
#' @param formula_full the model formula of the full model
#' @param formula_reduced the model formula of the model missing the covariate
#'  of interest
#' @param data data.frame with columnnames equal to variables specified in formulas
#' @param type regression type one of c("lm","glm","glmer")
#' @param family the family of the generalized models
#' @param weights the weights for the generalized models
#' @export
#' @details See \href{http://www.iaeng.org/publication/WCE2008/WCE2008_pp1063-1069.pdf}{
#'  A Signal-to-Noise Ratio Estimator for Generalized Linear Model Systems} for
#'  more info.
#' @examples
#'  glm_formula <- formula(Y ~ Surv(X,I) + Z)
#'  data <- simulate_data(100, glm_formula, type = "glm")
#'  glm_formula_reduced <- formula(Y ~ Z)
#'  glm_formula_full <- formula(Y ~ X + Z)
#'  general_signal_to_noise(
#'   glm_formula_full,
#'   glm_formula_reduced,
#'   data,
#'   type = "glm",
#'   weights = data$size_tot)
general_signal_to_noise <- function(formula_full, formula_reduced, data,
                                    type = c("lm","glm","glmer"),
                                    family = "binomial", weights){
  type <- match.arg(type)
  if (type == "lm"){
    args_full <- list(formula = formula_full, data = data)
    args_reduced <- list(formula = formula_reduced, data = data)
  } else if (type %in% c("glm","glmer")){
    args_full <- list(formula = formula_full, data = data, family = family, weights = weights)
    args_reduced <- list(formula = formula_reduced, data = data, family = family, weights = weights)
  }
  if (type == "glmer") type <- lme4::glmer
  fit_full <- do.call(type, args_full)
  fit_reduced <- do.call(type, args_reduced)
  snr <- (deviance(fit_reduced) - deviance(fit_full))/deviance(fit_full)
  return(snr)
}
