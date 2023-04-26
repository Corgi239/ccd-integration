CCD <- function(logden, 
                nParam, 
                logden.grad=NULL, 
                logden.hessian=NULL,
                mode=NULL, 
                th0=NULL,
                f0=1.1,  
                auto.scale=T, 
                optim.method='Nelder-Mead') {
  
  logden.count <- list("function"=0, "gradient"=0)
  
  ### Locate the mode of log-density and evaluate the Hessian
  if (is.null(mode)) {
    if (is.null(th0)) {
      th0 <- rep(0, nParam)
    }
    optim_res <- optim(th0, logden, gr = logden.grad, hessian = T, method=optim.method, control=list(fnscale=-1))
    mode <- optim_res$par
    H <- -optim_res$hessian
    logden.count <- optim_res$counts
  } else {
    if (is.null(logden.hessian)) {
      H <- -hessian(logden, mode)
    } else {
      H <- -logden.hessian(mode)
    }
  }
  
  ### Compute the whitening transformation from theta-space to z-space
  eig <- eigen(H)
  Lmbd <- diag(eig$values)
  V <- eig$vectors
  wt <- V %*% sqrtm(Lmbd)
  inv_wt <- solve(wt)
  
  ### Generate star points
  rad <- sqrt(nParam)
  star <- array(0, dim=c(nParam*2, nParam))
  for (i in 1:nParam) {
    star[2*i-1, i] <- rad
    star[2*i, i] <- -rad
  }
  
  ### Generate design points
  # Walsh indices to be included in design points (see Table II in Sanchez and Sanchez (2005))
  walsh <- c(1, 2, 4, 8, 15, 16, 32, 51, 64, 85, 106, 128, 150, 171, 219, 237,
             247, 256, 279, 297, 455, 512, 537, 557, 597, 643, 803, 863, 898,
             1024, 1051, 1070, 1112, 1169, 1333, 1345, 1620, 1866, 2048,
             2076, 2085, 2158, 2372, 2456, 2618, 2800, 2873, 3127, 3284,
             3483, 3557, 3763, 4096, 4125, 4135, 4176, 4435, 4459, 4469,
             4497, 4752, 5255, 5732, 5801, 5915, 6100, 6369, 6907, 7069,
             8192, 8263, 8351, 8422, 8458, 8571, 8750, 8858, 9124, 9314,
             9500, 10026, 10455, 10556, 11778, 11885, 11984, 13548, 14007,
             14514, 14965, 15125, 15554, 16384, 16457, 16517, 16609,
             16771, 16853, 17022, 17453, 17891, 18073, 18562, 18980,
             19030, 19932, 20075, 20745, 21544, 22633, 23200, 24167,
             25700, 26360, 26591, 26776, 28443, 28905, 29577, 32705)
  
  # determine the size of the full design matrix
  highest_index <- walsh[nParam] + 1
  folds <- ceiling(log2(highest_index))
  
  # construct the full factorial design matrix 
  H0 = 1
  for (i1 in 1:folds) {
    H0 <- rbind(cbind(H0, H0), cbind(H0, -H0))
  }
  
  # select points based on Walsh indices
  design <- H0[,1+walsh[1:nParam]]
  
  ### Join the star points, design points, and the mode
  z_points <- rbind(rep(0, nParam), star, design)
  point_count <- dim(z_points)[1]
  
  ### Shift the points along axial directions
  # matrix to store scaling coefficients for each coordinate of each point
  scaling_coeffs <- array(1, dim=dim(z_points))
  
  # compute log-density at the mode
  mode_density <- logden(mode)
  
  if (auto.scale) {
    # perform shifting for each main direction
    dir <- 1
    for (i in 1:(2*nParam)) {
      # Find the scaling parameter so that when we move sqrt(2) std
      # from the mode, the log density drops (approximately) by 1
      column_ind <- ceiling(i/2)
      probe_point <- rep(0, nParam)
      probe_point[column_ind] <- dir * sqrt(2)
      probe_density <- logden(mode + probe_point %*% inv_wt)
      logden.count[["function"]] <- logden.count[["function"]] + 1
      if (mode_density > probe_density) {
        scale <- sqrt(1/(mode_density - probe_density))
      } else {
        scale <- 1
      }
      
      # limit scaling to reasonable values
      scale <- pmax(pmin(scale, 10), 1/10)
      
      # scale coordinates of points in the given direction
      column <- z_points[, column_ind]
      scaling_coeffs[(column * dir) >= 0, column_ind] <- scale
      
      # flip the direction for next calculation
      dir <- -dir
    }
  }
  
  # apply the scaling coefficients
  z_points_scaled <- scaling_coeffs * z_points
  
  ### Return the points back to theta-space
  points <- f0 * z_points_scaled %*% inv_wt
  points <- points + mode[col(points)]
  
  ### Calculate CCD integration weights
  delta_i <- 1 / ((dim(points)[1] - 1) * (f0^2 - 1) * (1+exp(-nParam * f0^2 / 2)))
  delta_0 <- 1
  deltas <- append(c(delta_0), rep(delta_i, point_count - 1))
  
  ### Compute densities at integration points
  log.densities <- apply(points[-1,], 1, logden)
  log.densities <- append(mode_density, log.densities)
  log.densities <- log.densities - max(log.densities)
  densities <- exp(log.densities)
  
  ### Finalize the ingratiation weights
  weights <- deltas * densities
  weights <- weights / sum(weights)
  
  out <- list("points"=points, "weights"=weights, "densities"=densities, "logden.count"=logden.count)
  return(out)
}