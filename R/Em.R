
Em_tot_01_approx<- function(m, A, Pi, f0x, f1x,
                            x, eps, maxit, h) {
  diff <- list(eps +1, eps +1,eps +1,eps +1)
  i <- list(0)
  y <- as.list(rep(0,10))
  y[[6]] <- A
  y[[5]] <- list(gamma = matrix(Pi, nrow =1))
  y[[9]] <- f1x
  
  y[[10]] <- f0x
  while(((max(unlist(diff)) > eps) & (i[[1]] < maxit)))
  {
    
    y[[1]] <- y[[6]]
    # y[[2]] <- y[[5]]$gamma[1,]
    y[[3]] <- y[[10]]
    y[[4]] <- y[[9]]
    y[[5]] <- for_back(m, y[[6]], f0x=y[[3]], f1x = y[[4]], y[[5]]$gamma[1,])
    y[[6]] <- matrix(colSums(y[[5]]$ksi[-m,]) / rep(colSums(y[[5]]$gamma[-m,]),2),
                     byrow = FALSE, ncol = 2)
    y[[7]] <- density(x, weights = y[[5]]$gamma[,1]/sum(y[[5]]$gamma[,1]))
    y[[8]] <- density(x, weights = y[[5]]$gamma[,2]/sum(y[[5]]$gamma[,2]))
    y[[9]] <- approx(y[[8]]$x,y[[8]]$y,x)$y
    y[[10]] <- approx(y[[7]]$x,y[[7]]$y,x)$y
    
    diff[[1]] = max(abs( y[[6]] - y[[1]]))
    # diff[[2]] = max(abs(y[[5]]$gamma[1,] - y[[2]]))
    diff[[2]] = 0
    diff[[3]] = max(abs(y[[9]] - y[[4]]))
    diff[[4]] = max(abs(y[[10]] - y[[3]]))
    i[[1]] <- i[[1]]+1
    # gc()
  }
  
  b_f = for_back(m, y[[6]], y[[10]], y[[9]], y[[5]]$gamma[1,])
  return (list(A= y[[6]], Pi = y[[5]]$gamma[1,],
               fw_bc_EM= b_f,
               f1x =  y[[9]],
               f0x = y[[10]],
               i = i[[1]]))
}

Em_tot_approx<- function(m, A, Pi, f0x, f1x,
                         x, eps, maxit, h) {
  diff <- list(eps +1, eps +1,eps +1)
  i <- list(0)
  y <- as.list(rep(0,8))
  y[[6]] <- A
  y[[3]] <- Pi
  y[[8]] <- f1x
  while(((max(unlist(diff)) > eps) & (i[[1]] < maxit)))
  {
    # A_old = A
    # Pi_old = Pi
    # f1x_old = f1x
    # b_f = for_back(m, A, f0x, f1x, Pi)
    # sum_gamma = colSums(b_f$gamma[-m,])
    # sum_ksi = colSums(b_f$ksi[-m,])
    # Pi = b_f$gamma[1,]
    # A <- matrix(sum_ksi / rep(sum_gamma,2), byrow = FALSE, ncol = 2)
    # d0 <- density(x, weights = b_f$gamma[,1]/sum(b_f$gamma[,1]))
    # d1 <- density(x, weights = b_f$gamma[,2]/sum(b_f$gamma[,2]))
    # f1x <- approx(d1$x,d1$y,x)$y
    #
    # diff[1] = max(abs(A - A_old))
    # diff[2] = max(abs(Pi - Pi_old))
    # diff[3] = max(abs(f1x - f1x_old))
    
    y[[1]] <- y[[6]]
    # y[[2]] <-  y[[3]]
    y[[4]] <- y[[8]]
    y[[5]] <- for_back(m,y[[1]], f0x, y[[4]], y[[3]])
    y[[3]] <- y[[5]]$gamma[1,]
    y[[6]] <- matrix(colSums(y[[5]]$ksi[-m,]) / rep(colSums(y[[5]]$gamma[-m,]),2), byrow = FALSE, ncol = 2)
    y[[7]] <- density(x, weights = y[[5]]$gamma[,2]/sum(y[[5]]$gamma[,2]))
    y[[8]] <- approx(y[[7]]$x,y[[7]]$y,x)$y
    
    diff[[1]] = max(abs(y[[6]] - y[[1]]))
    # diff[[2]] = max(abs(y[[3]] - y[[2]]))
    diff[[2]] = 0
    diff[[3]] = max(abs(y[[8]] - y[[4]]))
    i[[1]] <- i[[1]]+1
    # gc()
    
  }
  
  b_f = for_back(m, y[[6]], f0x, y[[8]], y[[5]]$gamma[1,])
  return (list(A= y[[6]],Pi =  y[[5]]$gamma[1,],
               fw_bc_EM= b_f,
               f1x = y[[8]],
               f0x = f0x,
               i = i[[1]]))
}



Em <- function(m, A ,
               Pi,  f0x, f1x,
               x, eps = 0.0001,
               maxit =1000, h = 0.3, f0_known, approx = TRUE){
  if(approx){
    if(f0_known){
    a<-   Em_tot_approx(m, A ,
                    Pi,  f0x, f1x,
                    x, eps,
                    maxit, h)
    }else{
      Em_tot_01_approx(m, A ,
                       Pi,  f0x, f1x,
                       x, eps,
                       maxit, h)   
    }
  }else{
    if(f0_known){
      a<- Em_tot(m, A ,
             Pi,  f0x, f1x,
             x, eps,
             maxit, h)
      a$f0x <- f0x
      a
    }else{
      Em_tot_01(m, A ,
                Pi,  f0x, f1x,
                x, eps,
                maxit, h)   
    }  
  }
  
}