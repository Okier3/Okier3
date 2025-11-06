# Load libraries
library(geosphere)
library(dplyr)
setwd("D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/Eg/proj_01")
# Select BP and Control CSV files
bps_file <- file.choose()
control_file <- file.choose()

# Read CSV files
bps_data <- read.csv(bps_file, stringsAsFactors = FALSE)
control_data <- read.csv(control_file, stringsAsFactors = FALSE)

# Checking the CSVs for the required columns.
required_cols <- c("Name", "Latitude", "Longitude", "Orthometric_Height")

if (!all(required_cols %in% colnames(bps_data))) {
  stop("BPS CSV must contain columns: Name, Latitude, Longitude, Orthometric_Height")
}

if (!all(required_cols %in% colnames(control_data))) {
  stop("Control CSV must contain columns: Name, Latitude, Longitude, Orthometric_Height")
}

# Function to compute geodesic distance and azimuth using Vincenty's formula
compute_geodesic <- function(bps, control) {
  # Compute geodesic distance
  dist <- distVincentyEllipsoid(
    c(control$Longitude, control$Latitude),
    c(bps$Longitude, bps$Latitude)
  )
  
  # Compute initial azimuth (bearing) from control to BPS
  azimuth <- bearing(
    c(control$Longitude, control$Latitude),
    c(bps$Longitude, bps$Latitude)
  )
  # Convert negative azimuths to 0–360 range
  if (azimuth < 0) {
    azimuth <- azimuth + 360
  }
  
  # Compute ellipsoidal height difference
  height_diff <- bps$Orthometric_Height - control$Orthometric_Height
  
  return(data.frame(
    BPS_Name = bps$Name,
    BPS_Latitude = bps$Latitude,
    BPS_Longitude = bps$Longitude,
    BPS_Orthometric_Height = bps$Orthometric_Height,
    Control_Name = control$Name,
    Control_Latitude = control$Latitude,
    Control_Longitude = control$Longitude,
    Control_Orthometric_Height = control$Orthometric_Height,
    Distance_m = dist,
    Azimuth_deg = azimuth,
    Height_Diff_m = height_diff
  ))
}

# Compute distances and azimuths for all BP and control point pairs
results <- do.call(rbind, lapply(1:nrow(bps_data), function(i) {
  do.call(rbind, lapply(1:nrow(control_data), function(j) {
    compute_geodesic(bps_data[i, ], control_data[j, ])
  }))
}))

# Filter for controls within 50,000 meters
filtered_results <- results %>% filter(Distance_m < 50000)

# Display results
print(filtered_results)

# Save results to a CSV file
write.csv(filtered_results, "Filtered_Geodesic_Results.csv", row.names = FALSE)
#------------------------------------------------------------------------------
  #Height calculations.
#------------------------------------------------------------------------------
# Load required libraries
library(terra)  # Handling raster geoid models

# Load Geoid Models[5x5 arcminute resolution tidefree model]  
egm2008_file <- "D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/GGM_5/h_anomaly_EGM2008_5tf.GRD"
egm96_file <- "D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/GGM_5/h_anomaly_EGM96_5tf.GRD"
eigen_file <- "D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/GGM_5/h_anomaly_EIGEN-6C4_tf.GRD"
geco_file <- "D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/GGM_5/h_anomaly_GECO_tf.GRD"
sgg_file <- "D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/GGM_5/h_anomaly_SGG-UGM-1_tf.GRD"

geoid_egm2008 <- rast(egm2008_file)
geoid_egm96 <- rast(egm96_file)
geoid_eigen <- rast(eigen_file)
geoid_geco <- rast(geco_file)
geoid_sgg <- rast(sgg_file)
# Ensure CRS is WGS84 (EPSG:4210)
geoid_egm2008 <- project(geoid_egm2008, "EPSG:4210")
geoid_egm96 <- project(geoid_egm96, "EPSG:4210")
geoid_eigen <- project(geoid_eigen, "EPSG:4210")
geoid_geco <- project(geoid_geco, "EPSG:4210")
geoid_sgg <- project(geoid_sgg, "EPSG:4210")
# Function to Compute Geoid heights & Orthometric Heights
compute_heights <- function(df, geoid_egm2008, geoid_egm96,geoid_eigen,geoid_geco,geoid_sgg) {
  # Extract Geoid Heights from GRD
  geoid_height_egm2008 <- terra::extract(geoid_egm2008, df[, c("Longitude", "Latitude")], method = "bilinear")[, 2]
  geoid_height_egm96 <- terra::extract(geoid_egm96, df[, c("Longitude", "Latitude")], method = "bilinear")[, 2]
  geoid_height_eigen <- terra::extract(geoid_eigen, df[, c("Longitude", "Latitude")], method = "bilinear")[, 2]
  geoid_height_geco <- terra::extract(geoid_geco, df[, c("Longitude", "Latitude")], method = "bilinear")[, 2]
  geoid_height_sgg <- terra::extract(geoid_sgg, df[, c("Longitude", "Latitude")], method = "bilinear")[, 2]
  
  # Compute Ellipsoidal Heights(N=h-H)
  df$Geoid_Height_EGM2008 <- geoid_height_egm2008
  df$Geoid_Height_EGM96 <- geoid_height_egm96
  df$Geoid_Height_EIGEN <- geoid_height_eigen
  df$Geoid_Height_GECO <- geoid_height_geco
  df$Geoid_Height_SGG <- geoid_height_sgg
  df$Ellipsoidal_Height_EGM2008 <- df$Geoid_Height_EGM2008 +  df$Orthometric_Height
  df$Ellipsoidal_Height_EGM96 <-  df$Geoid_Height_EGM96 + df$Orthometric_Height
  df$Ellipsoidal_Height_EIGEN <-  df$Geoid_Height_EIGEN + df$Orthometric_Height
  df$Ellipsoidal_Height_GECO <-  df$Geoid_Height_GECO + df$Orthometric_Height
  df$Ellipsoidal_Height_SGG <-  df$Geoid_Height_SGG + df$Orthometric_Height
  
  return(df)
}

# Applying Function to BPs and Control Data
bps_data <- compute_heights(bps_data, geoid_egm2008, geoid_egm96,geoid_eigen,geoid_geco,geoid_sgg)
control_data <- compute_heights(control_data, geoid_egm2008, geoid_egm96,geoid_eigen,geoid_geco,geoid_sgg)

# Combine Results
final_results <- rbind(
  cbind(Type = "BPS", bps_data),
  cbind(Type = "Control", control_data)
)

# Results
print(final_results)

# Save to CSV
write.csv(final_results, "Geoid_Orthometric_Heights.csv", row.names = FALSE)
#------------------------------------------------------------------------------
#EGM2008
#------------------------------------------------------------------------------
library(terra)
library(dplyr)
library(sf)
library(MASS)
library(stats)

# Load First CSV Geoid orthometric Heights
geoid_file <- file.choose()  
geoid_data <- read.csv(geoid_file, stringsAsFactors = FALSE,header=TRUE)

# Load Second CSV filtered Geodesic file
gps_file <- file.choose()  
gps_data <- read.csv(gps_file, stringsAsFactors = FALSE, header=TRUE)
Name_BPS<-merge(geoid_data,gps_data, by.x="Name",by.y="BPS_Name", all.x=TRUE)

# Get Orthometric Heights from Geoid Data Without Merging
gps_data$BPS_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_EGM2008[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_EGM2008[match(gps_data$Control_Name, geoid_data$Name)]
gps_data$BPS_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$Control_Name, geoid_data$Name)]

# View the updated GPS data
head(gps_data)

# Compute Height Differences
gps_data <- gps_data %>%
  mutate(
    Delta_h = BPS_Ellipsoidal_Height - Control_Ellipsoidal_Height,  # Change in ellipsoidal height
    Delta_H = BPS_Orthometric_Height - Control_Orthometric_Height,  # Change in orthometric height
    S = Distance_m  # Horizontal distance
  )

# Compute Observation Vector R_h
gps_data$R_h <- ((gps_data$Delta_h - gps_data$Delta_H) / gps_data$S) * (-1)

# Convert Azimuths to Radians
gps_data$Azimuth_rad <- gps_data$Azimuth_deg * (pi / 180)

# Initialize results storage
results <- data.frame(BPS_Name = character(), Xi = numeric(), Eta = numeric(), 
                      Xi_arcsec = numeric(), Eta_arcsec = numeric(),
                      sigma0 = numeric(), sigma_xi = numeric(), sigma_eta = numeric(), stringsAsFactors = FALSE)

bps_list <- unique(gps_data$BPS_Name)

for (bps in bps_list) {
  # Subset data for the current BPS point
  subset_data <- gps_data %>% filter(BPS_Name == bps)
  # Ensure at least 2 control points for least squares solution
  if (nrow(subset_data) < 2) {
    cat("Skipping", bps, "- Not enough control points\n")
    next
  }
  print(subset_data)
  
  # Design matrix A
  A <- cbind(cos(subset_data$Azimuth_rad), sin(subset_data$Azimuth_rad))
  
  # Observation vector R
  R <- matrix(subset_data$R_h, ncol = 1)
  
  # Weight matrix W = diag(1/d^2)
  W <- diag(1 / (subset_data$S^2))
  
  # Weighted least squares solution: X = (Aᵀ W A)^(-1) Aᵀ W R
  AtW <- t(A) %*% W
  N <- AtW %*% A
  X <- solve(N) %*% (AtW %*% R)
  
  # Residuals
  v <- R - A %*% X
  dof <- nrow(A) - ncol(A)  # Degrees of freedom
  
  # A posteriori variance factor
  sigma0_sq <- as.numeric(t(v) %*% W %*% v) / dof
  sigma0 <- sqrt(sigma0_sq)
  
  # Covariance matrix of estimates
  cov_X <- sigma0_sq * solve(N)
  sigma_xi <- sqrt(cov_X[1, 1])
  sigma_eta <- sqrt(cov_X[2, 2])
  
  # Convert results to arcseconds
  xi <- X[1, 1]
  eta <- X[2, 1]
  xi_arcsec <- xi * (180/pi) * 3600
  eta_arcsec <- eta * (180/pi) * 3600
  
  # Store results
  results <- rbind(results, data.frame(
    BPS_Name = bps, Xi = xi, Eta = eta,
    Xi_arcsec = xi_arcsec, Eta_arcsec = eta_arcsec,
    sigma0 = sigma0,
    sigma_xi = sigma_xi * (180/pi) * 3600,
    sigma_eta = sigma_eta * (180/pi) * 3600
  ))
}
print(results)
write.csv(results, "DoV_Results_2008_weight.csv", row.names = FALSE)
#
#-------------------------------------------------------------------------------
#EGM96
#--------------------------------------------------------------------------------
# Load First CSV Geoid Heights
geoid_file <- file.choose()  
geoid_data <- read.csv(geoid_file, stringsAsFactors = FALSE,header=TRUE)

# Load Second CSV Geodesic file
gps_file <- file.choose()  
gps_data <- read.csv(gps_file, stringsAsFactors = FALSE, header=TRUE)
Name_BPS<-merge(geoid_data,gps_data, by.x="Name",by.y="BPS_Name", all.x=TRUE)

# Get Orthometric Heights from Geoid Data Without Merging
gps_data$BPS_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_EGM96[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_EGM96[match(gps_data$Control_Name, geoid_data$Name)]
gps_data$BPS_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$Control_Name, geoid_data$Name)]

# View the updated GPS data
head(gps_data)

# Compute Height Differences
gps_data <- gps_data %>%
  mutate(
    Delta_h = BPS_Ellipsoidal_Height - Control_Ellipsoidal_Height,  # Change in ellipsoidal height
    Delta_H = BPS_Orthometric_Height - Control_Orthometric_Height,  # Change in orthometric height
    S = Distance_m  # Horizontal distance
  )

# Compute Observation Vector R_h
gps_data$R_h <- ((gps_data$Delta_h - gps_data$Delta_H) / gps_data$S) * (-1)

# Convert Azimuths to Radians
gps_data$Azimuth_rad <- gps_data$Azimuth_deg * (pi / 180)

# Initialize results storage
results <- data.frame(BPS_Name = character(), Xi = numeric(), Eta = numeric(), 
                      Xi_arcsec = numeric(), Eta_arcsec = numeric(),
                      sigma0 = numeric(), sigma_xi = numeric(), sigma_eta = numeric(), stringsAsFactors = FALSE)

bps_list <- unique(gps_data$BPS_Name)

for (bps in bps_list) {
  # Subset data for the current BPS point
  subset_data <- gps_data %>% filter(BPS_Name == bps)
  # Ensure at least 2 control points for least squares solution
  if (nrow(subset_data) < 2) {
    cat("Skipping", bps, "- Not enough control points\n")
    next
  }
 print(subset_data)
 
  # Design matrix A
  A <- cbind(cos(subset_data$Azimuth_rad), sin(subset_data$Azimuth_rad))
  
  # Observation vector R
  R <- matrix(subset_data$R_h, ncol = 1)
  
  # Weight matrix W = diag(1/d^2)
  W <- diag(1 / (subset_data$S^2))
  
  # Weighted least squares solution: X = (Aᵀ W A)^(-1) Aᵀ W R
  AtW <- t(A) %*% W
  N <- AtW %*% A
  X <- solve(N) %*% (AtW %*% R)
  
  # Residuals
  v <- R - A %*% X
  dof <- nrow(A) - ncol(A)  # Degrees of freedom
  
  # A posteriori variance factor
  sigma0_sq <- as.numeric(t(v) %*% W %*% v) / dof
  sigma0 <- sqrt(sigma0_sq)
  
  # Covariance matrix of estimates
  cov_X <- sigma0_sq * solve(N)
  sigma_xi <- sqrt(cov_X[1, 1])
  sigma_eta <- sqrt(cov_X[2, 2])
  
  # Convert results to arcseconds
  xi <- X[1, 1]
  eta <- X[2, 1]
  xi_arcsec <- xi * (180/pi) * 3600
  eta_arcsec <- eta * (180/pi) * 3600
  
  # Store results
  results <- rbind(results, data.frame(
    BPS_Name = bps, Xi = xi, Eta = eta,
    Xi_arcsec = xi_arcsec, Eta_arcsec = eta_arcsec,
    sigma0 = sigma0,
    sigma_xi = sigma_xi * (180/pi) * 3600,
    sigma_eta = sigma_eta * (180/pi) * 3600
  ))
}

print(results)
write.csv(results, "DoV_Results_96_weight.csv", row.names = FALSE)
#--------------------------------------------------------------------------------
#EIGEN-6C4
#--------------------------------------------------------------------------------
# Load First CSV Geoid Heights
geoid_file <- file.choose()  
geoid_data <- read.csv(geoid_file, stringsAsFactors = FALSE,header=TRUE)

# Load Second CSV Geodesic file
gps_file <- file.choose()  
gps_data <- read.csv(gps_file, stringsAsFactors = FALSE, header=TRUE)
Name_BPS<-merge(geoid_data,gps_data, by.x="Name",by.y="BPS_Name", all.x=TRUE)

# Get Orthometric Heights from Geoid Data Without Merging
gps_data$BPS_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_EIGEN[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_EIGEN[match(gps_data$Control_Name, geoid_data$Name)]
gps_data$BPS_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$Control_Name, geoid_data$Name)]

# View the updated GPS data
head(gps_data)

# Compute Height Differences
gps_data <- gps_data %>%
  mutate(
    Delta_h = BPS_Ellipsoidal_Height - Control_Ellipsoidal_Height,  # Change in ellipsoidal height
    Delta_H = BPS_Orthometric_Height - Control_Orthometric_Height,  # Change in orthometric height
    S = Distance_m  # Horizontal distance
  )

# Compute Observation Vector R_h
gps_data$R_h <- ((gps_data$Delta_h - gps_data$Delta_H) / gps_data$S) * (-1)

# Convert Azimuths to Radians
gps_data$Azimuth_rad <- gps_data$Azimuth_deg * (pi / 180)

# Initialize results storage
results <- data.frame(BPS_Name = character(), Xi = numeric(), Eta = numeric(), 
                      Xi_arcsec = numeric(), Eta_arcsec = numeric(),
                      sigma0 = numeric(), sigma_xi = numeric(), sigma_eta = numeric(), stringsAsFactors = FALSE)

bps_list <- unique(gps_data$BPS_Name)

for (bps in bps_list) {
  # Subset data for the current BPS point
  subset_data <- gps_data %>% filter(BPS_Name == bps)
  # Ensure at least 2 control points for least squares solution
  if (nrow(subset_data) < 2) {
    cat("Skipping", bps, "- Not enough control points\n")
    next
  }
  print(subset_data)
  
  # Design matrix A
  A <- cbind(cos(subset_data$Azimuth_rad), sin(subset_data$Azimuth_rad))
  
  # Observation vector R
  R <- matrix(subset_data$R_h, ncol = 1)
  
  # Weight matrix W = diag(1/d^2)
  W <- diag(1 / (subset_data$S^2))
  
  # Weighted least squares solution: X = (Aᵀ W A)^(-1) Aᵀ W R
  AtW <- t(A) %*% W
  N <- AtW %*% A
  X <- solve(N) %*% (AtW %*% R)
  
  # Residuals
  v <- R - A %*% X
  dof <- nrow(A) - ncol(A)  # Degrees of freedom
  
  # A posteriori variance factor
  sigma0_sq <- as.numeric(t(v) %*% W %*% v) / dof
  sigma0 <- sqrt(sigma0_sq)
  
  # Covariance matrix of estimates
  cov_X <- sigma0_sq * solve(N)
  sigma_xi <- sqrt(cov_X[1, 1])
  sigma_eta <- sqrt(cov_X[2, 2])
  
  # Convert results to arcseconds
  xi <- X[1, 1]
  eta <- X[2, 1]
  xi_arcsec <- xi * (180/pi) * 3600
  eta_arcsec <- eta * (180/pi) * 3600
  
  # Store results
  results <- rbind(results, data.frame(
    BPS_Name = bps, Xi = xi, Eta = eta,
    Xi_arcsec = xi_arcsec, Eta_arcsec = eta_arcsec,
    sigma0 = sigma0,
    sigma_xi = sigma_xi * (180/pi) * 3600,
    sigma_eta = sigma_eta * (180/pi) * 3600
  ))
}
print(results)
write.csv(results, "DoV_Results_EIGEN_weight.csv", row.names = FALSE)
#--------------------------------------------------------------------------------
#GECO
#--------------------------------------------------------------------------------
# Load First CSV Geoid Heights
geoid_file <- file.choose()  
geoid_data <- read.csv(geoid_file, stringsAsFactors = FALSE,header=TRUE)

# Load Second CSV Geodesic file
gps_file <- file.choose()  
gps_data <- read.csv(gps_file, stringsAsFactors = FALSE, header=TRUE)
Name_BPS<-merge(geoid_data,gps_data, by.x="Name",by.y="BPS_Name", all.x=TRUE)

# Get Orthometric Heights from Geoid Data Without Merging
gps_data$BPS_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_GECO[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_GECO[match(gps_data$Control_Name, geoid_data$Name)]
gps_data$BPS_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$Control_Name, geoid_data$Name)]

# View the updated GPS data
head(gps_data)

# Compute Height Differences
gps_data <- gps_data %>%
  mutate(
    Delta_h = BPS_Ellipsoidal_Height - Control_Ellipsoidal_Height,  # Change in ellipsoidal height
    Delta_H = BPS_Orthometric_Height - Control_Orthometric_Height,  # Change in orthometric height
    S = Distance_m  # Horizontal distance
  )

# Compute Observation Vector R_h
gps_data$R_h <- ((gps_data$Delta_h - gps_data$Delta_H) / gps_data$S) * (-1)

# Convert Azimuths to Radians
gps_data$Azimuth_rad <- gps_data$Azimuth_deg * (pi / 180)

# Initialize results storage
results <- data.frame(BPS_Name = character(), Xi = numeric(), Eta = numeric(), 
                      Xi_arcsec = numeric(), Eta_arcsec = numeric(),
                      sigma0 = numeric(), sigma_xi = numeric(), sigma_eta = numeric(), stringsAsFactors = FALSE)

bps_list <- unique(gps_data$BPS_Name)

for (bps in bps_list) {
  # Subset data for the current BPS point
  subset_data <- gps_data %>% filter(BPS_Name == bps)
  # Ensure at least 2 control points for least squares solution
  if (nrow(subset_data) < 2) {
    cat("Skipping", bps, "- Not enough control points\n")
    next
  }
  print(subset_data)
  
  # Design matrix A
  A <- cbind(cos(subset_data$Azimuth_rad), sin(subset_data$Azimuth_rad))
  
  # Observation vector R
  R <- matrix(subset_data$R_h, ncol = 1)
  
  # Weight matrix W = diag(1/d^2)
  W <- diag(1 / (subset_data$S^2))
  
  # Weighted least squares solution: X = (Aᵀ W A)^(-1) Aᵀ W R
  AtW <- t(A) %*% W
  N <- AtW %*% A
  X <- solve(N) %*% (AtW %*% R)
  
  # Residuals
  v <- R - A %*% X
  dof <- nrow(A) - ncol(A)  # Degrees of freedom
  
  # A posteriori variance factor
  sigma0_sq <- as.numeric(t(v) %*% W %*% v) / dof
  sigma0 <- sqrt(sigma0_sq)
  
  # Covariance matrix of estimates
  cov_X <- sigma0_sq * solve(N)
  sigma_xi <- sqrt(cov_X[1, 1])
  sigma_eta <- sqrt(cov_X[2, 2])
  
  # Convert results to arcseconds
  xi <- X[1, 1]
  eta <- X[2, 1]
  xi_arcsec <- xi * (180/pi) * 3600
  eta_arcsec <- eta * (180/pi) * 3600
  
  # Store results
  results <- rbind(results, data.frame(
    BPS_Name = bps, Xi = xi, Eta = eta,
    Xi_arcsec = xi_arcsec, Eta_arcsec = eta_arcsec,
    sigma0 = sigma0,
    sigma_xi = sigma_xi * (180/pi) * 3600,
    sigma_eta = sigma_eta * (180/pi) * 3600
  ))
}
print(results)
write.csv(results, "DoV_Results_GECO_weight.csv", row.names = FALSE)
#--------------------------------------------------------------------------------
#SGG-UGM-1
#--------------------------------------------------------------------------------
# Load First CSV Geoid Heights
geoid_file <- file.choose()  
geoid_data <- read.csv(geoid_file, stringsAsFactors = FALSE,header=TRUE)

# Load Second CSV Geodesic file
gps_file <- file.choose()  
gps_data <- read.csv(gps_file, stringsAsFactors = FALSE, header=TRUE)
Name_BPS<-merge(geoid_data,gps_data, by.x="Name",by.y="BPS_Name", all.x=TRUE)

# Get Orthometric Heights from Geoid Data Without Merging
gps_data$BPS_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_SGG[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Ellipsoidal_Height <- geoid_data$Ellipsoidal_Height_SGG[match(gps_data$Control_Name, geoid_data$Name)]
gps_data$BPS_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$BPS_Name, geoid_data$Name)]
gps_data$Control_Orthometric_Height <- geoid_data$Orthometric_Height[match(gps_data$Control_Name, geoid_data$Name)]

# View the updated GPS data
head(gps_data)

# Compute Height Differences
gps_data <- gps_data %>%
  mutate(
    Delta_h = BPS_Ellipsoidal_Height - Control_Ellipsoidal_Height,  # Change in ellipsoidal height
    Delta_H = BPS_Orthometric_Height - Control_Orthometric_Height,  # Change in orthometric height
    S = Distance_m  # Horizontal distance
  )

# Compute Observation Vector R_h
gps_data$R_h <- ((gps_data$Delta_h - gps_data$Delta_H) / gps_data$S) * (-1)

# Convert Azimuths to Radians
gps_data$Azimuth_rad <- gps_data$Azimuth_deg * (pi / 180)

# Initialize results storage
results <- data.frame(BPS_Name = character(), Xi = numeric(), Eta = numeric(), 
                      Xi_arcsec = numeric(), Eta_arcsec = numeric(),
                      sigma0 = numeric(), sigma_xi = numeric(), sigma_eta = numeric(), stringsAsFactors = FALSE)

bps_list <- unique(gps_data$BPS_Name)

for (bps in bps_list) {
  # Subset data for the current BPS point
  subset_data <- gps_data %>% filter(BPS_Name == bps)
  # Ensure at least 2 control points for least squares solution
  if (nrow(subset_data) < 2) {
    cat("Skipping", bps, "- Not enough control points\n")
    next
  }
  print(subset_data)
  
  # Design matrix A
  A <- cbind(cos(subset_data$Azimuth_rad), sin(subset_data$Azimuth_rad))
  
  # Observation vector R
  R <- matrix(subset_data$R_h, ncol = 1)
  
  # Weight matrix W = diag(1/d^2)
  W <- diag(1 / (subset_data$S^2))
  
  # Weighted least squares solution: X = (Aᵀ W A)^(-1) Aᵀ W R
  AtW <- t(A) %*% W
  N <- AtW %*% A
  X <- solve(N) %*% (AtW %*% R)
  
  # Residuals
  v <- R - A %*% X
  dof <- nrow(A) - ncol(A)  # Degrees of freedom
  
  # A posteriori variance factor
  sigma0_sq <- as.numeric(t(v) %*% W %*% v) / dof
  sigma0 <- sqrt(sigma0_sq)
  
  # Covariance matrix of estimates
  cov_X <- sigma0_sq * solve(N)
  sigma_xi <- sqrt(cov_X[1, 1])
  sigma_eta <- sqrt(cov_X[2, 2])
  
  # Convert results to arcseconds
  xi <- X[1, 1]
  eta <- X[2, 1]
  xi_arcsec <- xi * (180/pi) * 3600
  eta_arcsec <- eta * (180/pi) * 3600
  
  # Store results
  results <- rbind(results, data.frame(
    BPS_Name = bps, Xi = xi, Eta = eta,
    Xi_arcsec = xi_arcsec, Eta_arcsec = eta_arcsec,
    sigma0 = sigma0,
    sigma_xi = sigma_xi * (180/pi) * 3600,
    sigma_eta = sigma_eta * (180/pi) * 3600
  ))
}
print(results)
write.csv(results, "DoV_Results_SGG_weight.csv", row.names = FALSE)
#-----------------------------------------------------------------------------------
#PLOTTING DISTRIBUTION
# Open PDF to save the plots
pdf("Network_Distribution_Plots.pdf", width = 7, height = 7)

bps_list <- unique(gps_data$BPS_Name)

for (bps in bps_list) {
  subset_data <- gps_data[gps_data$BPS_Name == bps, ]
  
  # Skip if no valid data
  if (nrow(subset_data) == 0) next
  
  # BPS coordinates (assuming same for all rows)
  bps_lat <- subset_data$BPS_Latitude[1]
  bps_lon <- subset_data$BPS_Longitude[1]
  
  # Create a combined data frame: control points and BPS
  controls <- data.frame(
    Name = subset_data$Control_Name,  # Assuming you have control point names
    Latitude = subset_data$Control_Latitude,
    Longitude = subset_data$Control_Longitude,
    Type = "Control Point"
  )
  
  bps_point <- data.frame(
    Name = bps,
    Latitude = bps_lat,
    Longitude = bps_lon,
    Type = "BPS"
  )
  
  network_points <- rbind(controls, bps_point)
  
  # Plot
  p <- ggplot() +
    geom_point(data = network_points, aes(x = Longitude, y = Latitude, shape = Type, color = Type), size = 3) +
    geom_text(data = network_points, aes(x = Longitude, y = Latitude, label = Name), hjust = -0.2, vjust = -0.5, size = 3) +
    geom_segment(data = controls, aes(x = bps_lon, y = bps_lat, xend = Longitude, yend = Latitude),
                 linetype = "dashed", color = "black") +
    scale_shape_manual(values = c("BPS" = 8, "Control Point" = 16)) +
    scale_color_manual(values = c("BPS" = "red", "Control Point" = "blue")) +
    labs(title = paste("Network Distribution for", bps),
         x = "Longitude", y = "Latitude") +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
}

dev.off()

