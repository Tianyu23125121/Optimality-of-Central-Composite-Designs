#################################################################################################
######### PROGRAM FOR SEARCHING EXACT OPTIMUM UNBLOCKED DESIGN FOR INFERENCE #################### 
#################################################################################################
source("functions_for_prediction_crit.R")

# ------------------------- PARAMETERS SETUP ---------------------------------------------------
# NUMBER OF FACTORS
K <- 6

# LEVELS OF EACH FACTOR                          
Levels <- list(-1:1, -1:1, -1:1, -1:1, -1:1,-1:1)

# DESIGN REGION: ENTER 'Y' FOR CUBIC OR 'N' FOR SPHERICAL     
Cubic <- 'N'

# INDICATORS OF TERMS IN THE SECOND ORDER MODEL
Terms <- matrix(c(rep(1, K), rep(1, K), rep(1, (K*(K-1)/2))), nr = 1)
Npar <- sum(Terms)

# NUMBER OF 'TRIES' FOR THE EXCHANGE ALGORITHM 
Ntries <- 50

# COMPOUND CRITERIA WEIGHTS (ID Criterion)
Kappa <- c(0, 0, 0, 0, 0, 0, 0, 1, 0)

# CONFIDENCE COEFFICIENTS
prob1 <- 0.95  
prob2 <- 0.95
prob3 <- 0.95

# BONFERRONI'S CORRECTION OF ALPHA FOR AP CRITERION?                                  
MC <- 'N'      
prob2 <- ifelse(MC == 'Y', prob2^(1/Npar), prob2)

# WEIGHTS FOR PARAMETERS (W DIAGONAL)
W <- matrix(rep(1, Npar), nr = 1)                              
for (i in (K+1):(2*K)) {
  if (Terms[i] == 1)
    W[sum(Terms[1:i])] <- W[sum(Terms[1:i])]/4
}
W <- matrix(W/sum(W), nr = 1)

##############################################################################################
###################### FUNCTION TO GENERATE CCD DESIGN #######################################
##############################################################################################
generateCCD <- function(K, center_points = 4, region = "spherical") {
  factorial_points <- expand.grid(rep(list(c(-1, 1)), K))
  factorial_points <- as.matrix(factorial_points[1:(2^(K-1)), ])
  
  alpha <- if (region == "spherical") sqrt(K) else 1.0
  axial_points <- do.call(rbind, lapply(1:K, function(i) {
    point_pos <- rep(0, K); point_neg <- rep(0, K)
    point_pos[i] <- alpha; point_neg[i] <- -alpha
    rbind(point_pos, point_neg)
  }))
  
  center_points_matrix <- matrix(rep(0, K * center_points), ncol = K, byrow = TRUE)
  ccd_design <- rbind(factorial_points, axial_points, center_points_matrix)
  return(ccd_design)
}

##############################################################################################
############################### ITERATE OVER DIFFERENT RUNS ##################################
##############################################################################################
# Set different run counts and number of center points
runs <- c(50,51,52,53,54,55)          # Total number of runs
center_points_list <- c( 6,7,8,9,10,11)  # Number of center points

results <- list()

for (i in 1:length(runs)) {
  cat("\n================== RUN:", runs[i], "CENTER POINTS:", center_points_list[i], "==================\n")
  
  # Set the current run count and number of center points
  N <- runs[i]
  center_points <- center_points_list[i]
  
  # 1. Run the exchange algorithm to get the optimized design
  design <- SearchTreat()
  opt_matrix <- matrixX(design$XBest, Terms)$X
  M_opt <- t(opt_matrix) %*% opt_matrix
  moment <- mat.I(Npar, K, Cubic)$S
  vol <- ifelse(Cubic == "Y", (2^K), ((pi^(K/2)) / gamma(K/2 + 1)) * sqrt(K)^K)
  opt_ID <- criterion(opt_matrix, M_opt, W, moment, Cubic, vol, Kappa)$critI_D
  
  # 2. Generate the CCD design and calculate the ID criterion
  ccd_design <- generateCCD(K, center_points, region = "spherical")
  ccd_matrix <- matrixX(ccd_design, Terms)$X
  M_ccd <- t(ccd_matrix) %*% ccd_matrix
  ccd_ID <- criterion(ccd_matrix, M_ccd, W, moment, Cubic, vol, Kappa)$critI_D
  
  # 3. Compare the results
  cat("  CCD Design ID Criterion       :", round(ccd_ID, 4), "\n")
  cat("  Optimized Design ID Criterion :", round(opt_ID, 4), "\n")
  
  if (opt_ID > ccd_ID) {
    cat("Result: The optimized design is better under ID criterion.\n")
  } else {
    cat("Result: The CCD design is better or equally optimal under ID criterion.\n")
  }
  
  # Save the current results
  results[[paste0("Run_", runs[i], "_Centers_", center_points)]] <- list(
    CCD_ID = ccd_ID,
    Optimized_ID = opt_ID,
    CCD_Design = ccd_design,
    Optimized_Design = design$XBest
  )
}

##############################################################################################
############################ PRINTING SUMMARY RESULTS ########################################
##############################################################################################
cat("\n\n================== SUMMARY OF RESULTS ==================\n")
for (name in names(results)) {
  cat("\n", name, "\n")
  cat("  CCD ID Criterion       :", round(results[[name]]$CCD_ID, 4), "\n")
  cat("  Optimized ID Criterion :", round(results[[name]]$Optimized_ID, 4), "\n")
  if (results[[name]]$Optimized_ID > results[[name]]$CCD_ID) {
    cat("  Result: Optimized design is better.\n")
  } else {
    cat("  Result: CCD design is better or equally optimal.\n")
  }
}
