

#' Title
#' 
#' @useDynLib hmm.sanssouci
#' @importFrom Rcpp evalCpp
#' @importFrom  magrittr "%>%"
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom tibble rowid_to_column
#' @importFrom tibble enframe
#' @importFrom tibble tibble
#' @param x numeric vector, of statistics (order using a given order)
#' @param al numeric, the risk
#' @param sel_function a function that return a tibble with selected set (see details)
#' @param delta anumeric,  risk delta to share the risk alpha between the bootstrap part and the estimated part
#' @param h numeric, the window size for the kde 
#' @param f0_known logical, wether f0 is known (if f0_known =TRUE the initialisation will be the true law under H0)
#' @param norm_init logical, wether the initialisation is normal or not (consider student)
#' @param n_boot numeric,  number of bootstrap sample
#' @param seuil numeric,  threshold for selected pvalues
#' @param m0_init numeric,  expectency under H0 (if norm_init = TRUE)
#' @param sd0_init numeric,  expectency under H0
#' @param df_init numeric, if norm_init = FALSE, student degree of freedom
#' @param max_pi0 the maximum value for the first estimation of pi0
#' @param approx wheter the kde is approximated using linear interpolation with a large range of values or calculate for every point.
#'
#' @return
#' @export
#'
#' @examples
#' m <- 2000
#' theta <- sim_markov(m, Pi = c(0.8,0.2), A = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T))
#' x <- rep(0, m)
#'  x[theta == 0] <- rnorm(sum(theta ==0))
#'  x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#'  sanssouci.hmm(x, al= 0.1, sel_function = Selection_delta) 
sanssouci.hmm <- function(x, al, sel_function= Selection_delta, delta = 0.1, h = 0.3,
                          f0_known = TRUE, 
                          norm_init = TRUE,
                          n_boot = 20,
                          seuil= 0.05,
                          m0_init = 0, sd0_init = 1, df_init = NULL,
                          max_pi0= 0.99999,
                          approx = TRUE, 
                          min_size =5, 
                          min_jump = 3, type_init = "given", drop_sel = TRUE){
  
Res <- res_all(x, al, sel_function, delta = 0.1, h = 0.3,
               f0_known = TRUE, 
               norm_init = TRUE,
               n_boot = 20,
               seuil= 0.05,
               m0_init = 0, sd0_init = 1, df_init = NULL,
               max_pi0= 0.99999,
               approx = TRUE, 
               min_size =5, 
               min_jump = 3, type_init = "given", drop_sel = TRUE)
Res$al <- al
Res$delta <- delta
return(mise_en_forme(Res))
  
} 


