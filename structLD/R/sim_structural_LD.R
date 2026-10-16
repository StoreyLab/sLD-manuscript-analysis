#' sim_structural_LD
#'
#' This function draws pairs of genotypes \eqn{X_a}, \eqn{X_b} in structural LD 
#' according to a multinomial distribution parameterized by individual-specific 
#' haplotype frequencies.
#'
#' @param pi_a Length-\eqn{n} vector of individual-specific allele frequencies (IAF) 
#'      at locus \eqn{a}. 
#' @param pi_b Length-\eqn{n} vector of IAFs at locus \eqn{b}. 
#' @param LD_sign Direction of LD (1 for positive LD, -1 for negative LD).
#' @param D_prime A value between 0 and 1 denoting the proportion of the maximum or
#'      minimum possible LD coefficient.
#'
#' @return Genotypes at loci \eqn{a} and \eqn{b} for \eqn{n} individuals. 
#'     
#' @note pi_a and pi_b must have the same length \eqn{n}. LD_sign must either equal 1 or -1. 
#'      D_prime must be in the interval \eqn{[0,1]}. 
#'
#' @examples
#' pi_a <- c(rep(0.25, 50), rep(0.3, 50))
#' pi_b <- c(rep(0.3, 50), rep(0.45, 50))
#' LD_sign <- 1
#' D_prime <- 0.4 
#' X_ab <- sim_structural_LD(pi_a, pi_b, LD_sign, D_prime)
#' 
#' @export

sim_structural_LD <- function(pi_a, pi_b, LD_sign, D_prime) {
  # Check arguments
  if (!is.numeric(pi_a) || !is.numeric(pi_b)) {
    stop("'pi_a' and 'pi_b' must be numeric vectors.")
  }
  if (length(pi_a) != length(pi_b)) {
    stop("'pi_a' and 'pi_b' must be of the same length.")
  }
  if (any(pi_a < 0 | pi_a > 1) || any(pi_b < 0 | pi_b > 1)) {
    stop("Allele frequencies must lie between 0 and 1.")
  }
  if (!LD_sign %in% c(1, -1)) {
    stop("'LD_sign' must equal either 1 or -1.")
  }
  if (!is.numeric(D_prime) || length(D_prime) != 1 || D_prime < 0 || D_prime > 1) {
    stop("'D_prime' must be a single numeric value between 0 and 1.")
  }
  
  # Compute the individual-specific LD coefficient for each individual 
  if (LD_sign == 1) {
    max_D <- pmin(pi_a * (1 - pi_b), (1 - pi_a) * pi_b)
    LD_offset <- D_prime * max_D 
  } else {
    min_D <- pmin(pi_a * pi_b, (1 - pi_a) * (1 - pi_b))
    LD_offset <- -D_prime * min_D
  }
  
  # Compute individual-specific haplotype frequencies
  pi_ab <- pi_a*pi_b + LD_offset
  pi_aB <- pi_a*(1 - pi_b) - LD_offset
  pi_Ab <- (1 - pi_a)*pi_b - LD_offset
  pi_AB <- (1 - pi_a)*(1 - pi_b) + LD_offset
  haplo_freqs <- cbind(pi_ab, pi_aB, pi_Ab, pi_AB)
  
  # Draw haplotypes and compute genotypes 
  n_ind <- length(pi_a)
  haplo <- t(sapply(1:n_ind, function(i) rmultinom(1, size = 2, prob = haplo_freqs[i, ])))
  X_a <- haplo[, 1] + haplo[, 2]
  X_b <- haplo[, 1] + haplo[, 3]
  X_ab <- rbind(X_a, X_b)
  
  return(X_ab)
}