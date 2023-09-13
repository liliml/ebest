#title in line below
#' Empirical Bayes Estimators for Population-Specific Association
#'
#description in line below
#' ebest provides joint and pairwise empirical Bayes (EB) estimators for population-specific association analysis,
#' while leveraging information from other populations, accounting for both the uncertainty of the estimates and
#' heterogeneity across populations. ebest also provides the conventional meta-analysis fixed effect estimator
#' for combining the estimates across populations.
#'
#' @param K Number of populations
#' @param betas Vector of regression coefficient estimates for the K populations, the first one is considered as the population of interest
#' @param ses Vector of standard errors of the corresponding regression coefficient estimates for the \emph{K} populations
#'
#' @details Estimates and standard errors for the meta-analysis estimator (b.meta, se.meta), joint EB estimator (b.joint, se.joint),
#' and pairwise EB estimator (b.pairwise, se.pairwise).
#' If (K = 2), the joint and pairwise EB estimators are same and only the joint EB estimator is outputed.
#' @references Hsu, L., Kooperberg, A., Reiner, A.P. and Kooperberg, C., (2023) An empirical Bayes approach to improving population‐specific genetic association estimation by leveraging cross‐population data. Genetic Epidemiology, 47(1), pp.45-60.
#' @author Contributors: Li Hsu, Lilian Law
#' @note In the example give below, 3 populations are used (K value), and the 3 values for the population sizes, p.freq, and betas each given in the 3 vectors are used to calculate the parameters "beta" and "ses" in the function ebest.
#'
#' @examples
#'
#' K = 3
#' est = c(0.1, 0.2, 0.3)
#' se = c(0.2, 0.2, 0.2)
#' ebest(K = K, betas = est, ses = se)
#'
#' @import quadprog
#' @import tibble
#' @import stats
#' @export
ebest <- function(K, betas, ses) {
##### The program yielded three estimates: meta-analysis fixed effects; joint EB; pairwise EB
##### Input
#####    K: number of populations
#####    betas: vector of regression coefficient estimates for the K populations, the first one is considered as the population of interest
#####    ses: vector of standard errors of the corresponding regression coefficient estimates for the K populations
##### Output
#####    estimates and standard errors for the meta-analysis estimator (b.meta, se.meta)
#####    joint EB estimator (b.joint, se.joint)
#####    pairwise EB estimator (b.pairwise, se.pairwise)
#####    if (K=2), the joint and pairwise EB estimators are same and only the joint EB estimator is outputed

  if (K<2|length(betas)!=K|length(ses)!=K) {
    print('Check the input')
    output <- NA
    return(output)
  }
### meta estimator
  denom = sum(1/ses^2)
  meta.weights = 1/ses^2 / denom
  b.meta = sum(meta.weights * betas)
  se.meta = denom^(-0.5)

### joint eb estimator
  v = sum(meta.weights^2 * (betas - betas[1])^2)
  wt.joint = v / (v + ses[1]^2)
  b.joint = wt.joint*betas[1] + (1-wt.joint)*b.meta

  c1 = 1 - wt.joint
  c2 = 2*ses[1]^2*(b.meta-betas[1]) / (v + ses[1]^2)^2
  c3 = sum(meta.weights^2 * (betas - betas[1]))

  # added square
  se.joint = (ses[1]^2 * (
    (1 + (meta.weights[1]-1)*c1) + c2*c3)^2 +
      sum(ses[-1]^2 * (meta.weights[-1] * c1 - c2*meta.weights[-1]^2 *
                         (betas[-1] - betas[1]))^2))^(0.5)

### pairwise eb estimators

  if (K>2){
    meta.weights = ses[-1]^2/(ses[-1]^2 + ses[1]^2)
    b.metas = meta.weights*betas[1] + (1-meta.weights)*betas[-1]   ## pairwise meta-estimator
    ses.meta = ses[1]*ses[-1]/(ses[1]^2 + ses[-1]^2)^0.5

    eb.weights = (betas[1]-b.metas)^2 / ((betas[1]-b.metas)^2 + ses[1]^2)
    b.ebs = eb.weights*betas[1] + (1-eb.weights)*b.metas
    mu = betas[1] - b.metas
    d1 = ses[1]^4/(ses[1]^2 + ses[-1]^2)
    d2 = (ses[1]^2 - mu^2)/(ses[1]^2 + mu^2)^2
    se.ebs =
      (ses[1]^2 + (d1)^2*d2^2*(ses[1]^2+ses[-1]^2) - 2*d1*d2*ses[1]^2)^0.5
    # covariance matrix
    cov = cbind(1 - d1*d2) %*% rbind(1 - d1*d2) * ses[1]^2
    diag(cov) <- se.ebs^2

    cor.eb = cov2cor(cov)

    # minimize variance
    qp = quadprog::solve.QP(cov,
                  array(0, dim = c(K-1,1)),
                  cbind(rep(1, K-1), diag(K-1)),
                  array(c(1, rep(0, K-1)), dim = c(K,1)))
    wts = qp$solution
    se.pairwise = (2*qp$value)^0.5
    b.pairwise = sum(wts * b.ebs)
  }

  if (K==2){
    output <- tibble::tibble(
      b.meta = b.meta,
      se.meta = se.meta,
      b.joint = b.joint,
      se.joint = se.joint)
  }
  if (K > 2){
    output <- tibble::tibble(
      b.meta = b.meta,
      se.meta = se.meta,
      b.joint = b.joint,
      se.joint = se.joint,
      b.pairwise = b.pairwise,
      se.pairwise = se.pairwise)
  }

  return(output)
}

