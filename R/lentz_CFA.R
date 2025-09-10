lentz_CFA <- function(x){
  
  eps <- 1e-12
  tiny <- 1e-30
  max_iter <- 10000
  
  # Initialization:
  f <- b0 <- x
  C <- f
  if (abs(C) < tiny) C <- tiny
  D <- 0
  delta <- 0
  j <- 1
  
  repeat {
    # Continued fraction coefficients:
    if(j %% 2 == 1) {
      a <- ((j - 1) / 2)^2
      b <- 1
    } else {
      a <- - (j / 2)^2
      b <- x
    }
    
    # Lentz's algorithm update:
    D <- b + a * D
    if(abs(D) < tiny) D <- tiny
    D <- 1 / D
    
    C <- b + a / C
    if(abs(C) < tiny) C <- tiny
    
    delta <- C * D
    f <- f * delta
    
    # Check convergence:
    if(abs(delta - 1) < eps || j >= max_iter) break
    
    j <- j + 1
  }
  
  ans <- 1/f
  
  return(ans)
  
}