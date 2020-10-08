#' Calculate scales for specific censoring rates and covariance dependent censoring
#'
#' @param censoring_val censoring rates to obtain
#' @param log_ratio_val log ratios of censoring of two distributions. e.g.
#'  log(20/10) means 20 censored values because of covariate one and 10 censored
#'  values because of covariate two.
#' @param seq_to_test parameter space
#' @param n number of random samples for each parameter combination
#' @param shapes the shapes of the weibull distributions
#' @param scale_x scale of the true weibull distribution
#'
#' @return tibble
#' @export
#'
#' @examples scales_for_censoring(c(0.3,0.5),c(0,0.7))
scales_for_censoring <- function(
  censoring_val = c(0.3,0.5),
  log_ratio_val = c(0,0.7),
  seq_to_test = c(seq(0.01,0.99,length.out = 50), seq(1,10,by=0.5)),
  n = 10000,
  shapes = list(shape_x = 0.5,
                 shape_c1 = 1,
                 shape_c2 = 0.5),
  scale_x = 0.25){

  # without covariate dependent censoring (faster)
  if (all(log_ratio_val==0)){
    tmpgrid_scale <- expand.grid(seq_to_test)
    tmpgrid_ind <- expand.grid(seq_along(seq_to_test))
    outmat_cens <- matrix(NA,ncol=length(seq_to_test),nrow = length(seq_to_test))
    rownames(outmat_cens) <- round(seq_to_test,2)
    colnames(outmat_cens) <- round(seq_to_test,2)
    # integrate this function to get expected censoring rate, basically the probability
    # that one random value is smaller than another one
    f <- function(x,shape_x=10,scale_x=4,shape_c=1,scale_c=1){
      pweibull(x,shape_c,scale_c)*dweibull(x,shape_x,scale_x)
    }
    outmat_cens <- sapply(seq_len(length(seq_to_test)), function(x){
      integrate(f,
                shape_x=shapes[["shape_x"]],scale_x=scale_x,
                shape_c=shapes[["shape_c1"]],scale_c=tmpgrid_scale[x,1],
                lower=0,upper=Inf)$value
    })
    names(outmat_cens) <- round(seq_to_test,2)
    cond_to_test <- expand.grid(censoring_val = censoring_val)
    scales <- purrr::pmap(cond_to_test, function(censoring_val){
      run <- TRUE
      tol <- 0.01
      while (run){
        wind_new <- which(((outmat_cens>(censoring_val-tol)) &
                         (outmat_cens < (censoring_val+tol))),arr.ind = TRUE)
        if (length(wind_new) == 0){
          run <- FALSE
        } else {
          outmat_cens <- outmat_cens[wind_new]
          wind <- wind_new
          tol <- tol*0.9
        }
      }

      if(!exists("wind")) stop("No valid parameters found, increase parameter space",call. = FALSE)
      tibble::tibble(censoring = censoring_val, log_ratio = log_ratio_val,
                     C1 = as.double(names(outmat_cens)[wind[1]]), C2 = C1)
    })

    # with covariate dependent censoring (slow)
  } else {
    tmpgrid_scale <- expand.grid(seq_to_test,seq_to_test)
    tmpgrid_ind <- expand.grid(seq_along(seq_to_test),seq_along(seq_to_test))
    outmat_logratio <- matrix(NA,ncol=length(seq_to_test),nrow = length(seq_to_test))
    outmat_cens <- matrix(NA,ncol=length(seq_to_test),nrow = length(seq_to_test))
    rownames(outmat_logratio) <- round(seq_to_test,2)
    colnames(outmat_logratio) <- round(seq_to_test,2)
    rownames(outmat_cens) <- round(seq_to_test,2)
    colnames(outmat_cens) <- round(seq_to_test,2)

    purrr::walk(seq_len(length(seq_to_test)*length(seq_to_test)), function(x){
      xw <- rweibull(n,shapes[["shape_x"]],scale_x)
      c1w <- rweibull(n,shapes[["shape_c1"]],tmpgrid_scale[x,1])
      c2w <- rweibull(n,shapes[["shape_c2"]],tmpgrid_scale[x,2])
      c1c <- xw > c1w
      c2c <- xw > c2w
      c <- c(c1w[seq(1,round(n/2))],c2w[seq(round(n/2)+1,n)])
      xind <- xw > c
      log_ratio <- log(sum(c1c[seq(1,round(n/2))])/sum(c2c[seq(round(n/2)+1,n)]))
      outmat_logratio[tmpgrid_ind[x,1],tmpgrid_ind[x,2]] <<- log_ratio
      outmat_cens[tmpgrid_ind[x,1],tmpgrid_ind[x,2]] <<- sum(xind)/n
    })
    tolerance_log_ratio <- log(1.2)
    cond_to_test <- expand.grid(censoring_val = censoring_val,log_ratio_val = log_ratio_val)
    scales <- purrr::pmap(cond_to_test, function(censoring_val,log_ratio_val){
      wind <- which(((outmat_cens>(censoring_val-0.01)) &
                       (outmat_cens < (censoring_val+0.01)) &
                       (abs(outmat_logratio) > log_ratio_val) &
                       (abs(outmat_logratio) < (log_ratio_val+tolerance_log_ratio))),arr.ind = TRUE)
      if(dim(wind)[1] == 0) stop("No valid parameters found, increase parameter space",call. = FALSE)
      wind[ ,1] <- as.double(rownames(outmat_cens)[wind[,1]])
      wind[ ,2] <- as.double(rownames(outmat_cens)[wind[,2]])
      rownames(wind) <- NULL
      tibble::tibble(censoring = censoring_val, log_ratio = log_ratio_val, C1 = wind[1,1], C2 = wind[1,2])
    })
  }
  scales <- dplyr::bind_rows(scales)

  return(scales)
}
