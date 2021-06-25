# 
borne <- function( type_borne, sel, a, alpha, m, 
                   pval,
                   C, ZL, leaf_list){
  l_sel <- length(sel)
  if(l_sel== 0){ return(NA)}
  if(type_borne == "HMM"){
    return( min(which(a >= 1 - alpha) -1, l_sel))
  }
  
  if(type_borne == "HMM_small"){
    
    return( max(c(which(a < alpha),1)) - 1)
  }
  if(type_borne == "Simes"){
    # pS <- posthocBySimes(pval, sel, alpha, Rcpp = FALSE, verbose = FALSE) 
    #
    return(min(sapply(1:l_sel, function(k){ sum(pval[sel] > alpha * k / m) +
        k - 1})))
    
    # return(length(sel)- pS)
  }
  if(type_borne=="DKW_tree"){
    return( V.star(sel, C, ZL, leaf_list))
    
    
  }
  if(type_borne=="DKW"){
    # return( V.star(sel, C, ZL, leaf_list))
    C <- sqrt(0.5 * log(1 / alpha))
    pval_sel <- pval[sel]
    pi <- pval_sel[order(pval_sel)]
    DKW_fun <- function(i) {
      (C / (2 * (1 - pi[i])) +
         (C^2 / (4 * (1 - pi[i])^2) +
            (l_sel - i) / (1 - pi[i]))^(1 / 2))^2
    }
    return(min(c(sapply(1:l_sel, DKW_fun), l_sel)))
    
  }
}




## function to get different bounds for a vector of probabilities 
## the function getprobs calculate the quantiles for different probabilities using the matrix B.
get_probs <- function(sel, li0, Pis, f0x, f1x, probs){
  if(length(sel)< 2) {
    return(NA)
  }else{
    Pis_sel <- lapply(1:(length(sel)-1), function(i){
      Reduce("%*%", Pis[sel[i]: (sel[i +1]-1)])
    })
    l_sel <- length(sel)
    probs_ord <- sort(probs)
    petit_grand <- rep(0, length(probs))
    petit_grand[probs_ord <0.5] <- 1
    Probs <- getprobs( length(sel),  
                       probs_ord, 
                       length(probs),
                       petit_grand,
                       li0 = li0[sel], 
                       f0x[sel],
                       f1x[sel],
                       Pis_sel)
    Probs_ord <- Probs[order(order(probs))]
    return(Probs_ord)
  }
}