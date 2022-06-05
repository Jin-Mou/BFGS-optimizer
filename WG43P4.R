#### Work Group 43 by Yixin Wang and Jin Mou
#### Github repo: https://github.com/cxszry/gas3.git

#### CODE OVERVIEW
### function "bfgs" is built to implement the BFGS quasi-Newton 
### minimization method. in the function:
### (1) we set and calculate initial value.
### (2) we calculate the initial gradient using initial parameters.
### (3) we use while loop to search for the optimal objective.
### (4) we calculate Hessian matrix H by finite differencing.


#### define bfgs
bfgs <- function(theta, f, ..., tol=1e-5, fscale=1, maxit=100){
  ### This function is used to implement the BFGS quasi-Newton minimization method
  ### Input parameters explanation:
  ## theta: a vector of initial values for the optimization parameters
  ## f: the objective function to minimize:
  ##    its first argument is the vector of optimization parameters.
  ##    its second argument is a logical parameter. TRUE if gradients need to be computed
  ## ...: any arguments of f after the first two
  ## tol: the convergence tolerance
  ## fscale: a rough estimate of the magnitude of f at the optimum
  ## maxit: the maximum number of BFGS iterations to try
  
  ### (1) set and calculate initial value
  # build function fx which inherits parameters from function f
  fx <- function(theta) f(theta, ...) 
  g <- array(0, length(theta)) # create grad to store gradients
  # set the initial inverse Hessian matrix B to identity matrix
  B <- diag(length(theta)) 
  par_num <- length(theta) # the number of parameters
  interval <- 1e-7 # finite difference interval
  iter <- 0 # the initial number of iteration
  
  ### (2) get the initial gradient
  ## if function f does not give fx an attribute "gradient",
  ## then we calculate it by finite differencing, else we inherit it
  # no_grad is TRUE if gradient needs to be computed
  no_grad <- is.null(attr(fx(theta), "gradient"))
  f0 <- fx(theta) # operate the function fx using initial parameters
  if (no_grad) {
    for (i in 1:par_num) { # loop over parameters
      theta1 <- theta
      theta1[i] <- theta1[i] + interval # increase theta[i] by interval
      f1 <- fx(theta1) # compute f1 using parameters in theta1
      g[i] <- (f1 - f0)/interval # approximate the gradient
    }
  } else { # inherit gradients from fx
    g <- attr(fx(theta), "gradient")
  }
  g <- matrix(g, par_num, 1) # change the grad from array to matrix
  
  ### (3) while loop to search for the optimal objective:
  ### if the gradient vector is close enough to 0, which means we roughly reach
  ### the optimal objective, then break the loop.
  while (max(abs(g)) >= (abs(f0)+fscale)*tol) {
    ## (3.1) get the Quasi-Newton step that satisfies the second Wolfe condition
    delta <- -B %*% g # calculate the descent direction delta
    step_length = 1 # initial step length
    # using repeat to adjust step_length to find the appropriate step length.
    # if two following conditions are satisfied, then exit the loop
    repeat {
      g1 <- g # the next potential gradient
      theta2 <- theta + step_length*t(delta) # the new theta
      f2 <- fx(theta2) # compute f2 using parameters in theta2
      # calculate next gradient g1
      if (no_grad) {
        for (i in 1:par_num) { # loop over parameters
          theta3 <- theta2
          theta3[i] <- theta3[i] + interval # increase theta2[i] by interval
          f3 <- fx(theta3) # compute f3 using parameters in theta3
          g1[i] <- (f3 - f2)/interval # approximate the new gradient
        }
      } else { # inherit gradients from fx
        g1 <- attr(fx(theta2), "gradient")
        g1 <- matrix(g1, par_num, 1) # change the gradient from array to matrix
      }
      
      # if two conditions as follows are satisfied, then exit loop
      conditions_satisfied <- 0 # count the number of conditions satisfied
      # condition 1
      # increase the step length if the second Wolfe condition is not satisfied
      if (t(g1) %*% delta < 0.9*t(g) %*% delta) {
        step_length = 1.8 * step_length
        next
      } else {
        conditions_satisfied <- conditions_satisfied + 1
      }
      # condition 2
      # reduce the step length if it leads to bigger objective value
      if (f2 > f0) {
        step_length = 0.2 * step_length
        next
      } else {
        conditions_satisfied <- conditions_satisfied + 1
      }
      # if both conditions are satisfied, then break 
      if (conditions_satisfied == 2) break
    }

    ## (3.2) update parameters
    # update inverse Hessian matrix B using BFGS method
    s = matrix(theta2 - theta) # transpose of theta2-theta
    y = g1 - g
    m <- 1 / (t(s) %*% y)
    m <- as.numeric(m) # m should be calculated as a number
    B <- (diag(par_num) - m*s%*%t(y)) %*% B %*% (diag(par_num) - m*y%*%t(s)) +
      m*s%*%t(s)
    # update theta
    theta <- theta2
    # update gradient g
    g <- g1
    # update f0
    f0 <- fx(theta)
    # update iter
    iter <- iter + 1
  }
  
  ### (4) calculate Hessian matrix H
  H <- matrix(0,par_num,par_num) # finite diference Hessian
  if (no_grad) {
    for (i in 1:par_num) { # loop over parameters
      theta1 <- theta
      theta1[i] <- theta1[i] + interval # increase theta1[i] by interval
      ## calculate gradient of theta1
      f1 <- fx(theta1)
      for (j in 1:par_num) {
        theta2 <- theta1
        theta2[j] <- theta2[j] + interval # increase theta2[i] by interval
        f2 <- fx(theta2)
        g1[j] <- (f2 - f1)/interval # approximate the gradient
      }
    H[i, ] <- (g1 - g)/interval # approximate the Hessian matrix H
    }
  } else { # inherit gradients from fx
    for (i in 1:par_num) { # loop over parameters
      theta1 <- theta
      theta1[i] <- theta1[i] + interval # increase theta1[i] by interval
      ## calculate gradient of theta1
      g1 <- attr(fx(theta1), "gradient")
      g1 <- matrix(g1, par_num, 1) # change the gradient from array to matrix
      H[i, ] <- (g1 - g)/interval # approximate the Hessian matrix H
    }
    H <- 0.5 * (t(H) + H) # fix H to be symmetric
  }
  
  list(f=f0[1], theta=array(theta), g=array(g), iter=iter, H=H)
}  
