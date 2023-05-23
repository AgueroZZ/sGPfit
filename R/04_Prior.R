#' Compute the correction factor based on the PSD.
#'
#' This function outputs the correction factor that converts the prior on
#' PSD to the prior on the original SD, based on the prediction unit d and
#' frequency parameter a.
#'
#' @param a A positive scalar represents the frequency parameter.
#' @param d A positive scalar represents the prediction unit.
#' @return A scalar denotes the correction factor to recover the original SD.
#' @export
compute_d_step_sGPsd <- function(d,a){
  sqrt((1/(a^2))*((d/2) - (sin(2*a*d)/(4*a))))
}












