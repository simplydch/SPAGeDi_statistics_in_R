# SPAGeDi functions R

# Function to calculate distance between pairs using Spherical law of cosines which
# matches SpageDi, even though it says it uses Euclidian distances
# The Haversine formula produces similar results, the distm function from the geosphere
# can also produce either of these without transforming the coordinates

calc_spatial_dist <- function(data_file){
  pairs <- t(combn(data_file[,"id"], 2))
  ind_i_coords <- data_file[match(pairs[, 1],  data_file[,"id"]),c("longitude","latitude")] * pi/180
  ind_j_coords <- data_file[match(pairs[, 2],  data_file[,"id"]),c("longitude","latitude")] * pi/180
  
  # Sets locations to 0 when they are the same, this solves a issues where acos can return nan
  # When coordinates are the same due to floating point errors
  same_loc <- abs(ind_i_coords[,1]-ind_j_coords[,1]) < 1e-10 & 
    abs(ind_i_coords[,2]-ind_j_coords[,2]) < 1e-10
  ind_i_coords[same_loc & !is.na(same_loc),] <- 0
  ind_j_coords[same_loc & !is.na(same_loc),] <- 0
  
  #ifelse(ind_i_coords[,1]==ind_j_coords[,1] & ind_i_coords[,2]==ind_j_coords[,2], 0,
  spatial_dists <- acos(sin(ind_i_coords[,2])*sin(ind_j_coords[,2])+cos(ind_i_coords[,2])*cos(ind_j_coords[,2])*
                          cos(ind_j_coords[,1]-ind_i_coords[,1]))*6371
  
  return(data.frame(id1 = pairs[,1], id2=pairs[,2], dist=spatial_dists))
}

# Function to calculate Moran's i


calc_morans_i <- function(IDs, genotypes, spatial_distances) {
  
  genotypes <- genotypes/2
  num_non_na_per_SNP <- colSums(genotypes * 0 + 1, na.rm = TRUE)
  
  
  allele_freq <- colSums(genotypes, na.rm = TRUE) / num_non_na_per_SNP
  
  # Calculate per SNP population variance
  
  vars <- rowSums((t(genotypes) - colMeans(genotypes, na.rm = T))^2, na.rm = T) / num_non_na_per_SNP
  
  
  valid <- num_non_na_per_SNP - 1
  
  
  pairs <- t(combn(IDs, 2))
  bias_cor <- vars / valid
  
  ind_i_ref <- match(pairs[, 1], IDs)
  ind_j_ref <- match(pairs[, 2], IDs)
  
  
  out_moranI <- sapply(1:nrow(pairs), function(x) {
    ind_i <- ind_i_ref[x]
    ind_j <- ind_j_ref[x]
    
    geno <- genotypes[c(ind_i, ind_j), ]
    y <- complete.cases(t(geno))
    
    geno_complete <- geno[, y]
    
    allele1 <- (geno_complete[1, ] - allele_freq[y]) * (geno_complete[2, ] - allele_freq[y])
    sum(allele1 + (bias_cor[y])) / ((sum(vars[y])))
  })
  
  return(out_moranI)
}


# Function to calculate correlograms

create_correlogram <- function(dist_classes, spatial_distances, MoransI_values) {
  spatial_distances[spatial_distances == -1] <- 0
  dist_mat <- cbind(MoransI_values, spatial_distances)
  
  correlogram <- sapply(2:length(dist_classes), function(x) {
    low <- dist_classes[x - 1]
    high <- dist_classes[x]
    mat <- dist_mat[dist_mat[, 2] > low & dist_mat[, 2] <= high, , drop = FALSE]
    if (!is.null(mat) && nrow(mat) > 0) {
      mean(mat[, 1], na.rm=TRUE)
    } else {
      NA
    }
  })
  return(correlogram)
}

# Function to calculate  Roussets a

calc_roussets_a <- function(IDs, genotypes) {
  
  genotypes <- genotypes/2
  
  pairs <- t(combn(IDs, 2))
  
  ind_i_ref <- match(pairs[, 1], IDs)
  ind_j_ref <- match(pairs[, 2], IDs)
  
  temp_geno <- genotypes * 0
  temp_geno[genotypes==0.5] <- 1
  
  num_non_na_per_loc <- colSums(genotypes * 0 + 1, na.rm = TRUE)
  geno_sums <- colSums(temp_geno,na.rm=T)
  rm(temp_geno)
  
  
  out_rousset_a <- sapply(1:nrow(pairs), function(x) {
    ind_i <- ind_i_ref[x]
    ind_j <- ind_j_ref[x]
    
    
    geno <- genotypes[c(ind_i, ind_j), ]
    y <- complete.cases(t(geno))
    geno_complete <- geno[, y]
    
    Hwt <- sum(geno_complete==0.5)/(ncol(geno_complete)*2)
    Hot <- mean(geno_sums[y]/num_non_na_per_loc[y])
    #Hbt <- sum(abs(geno_complete[1,] - geno_complete[2,]))/ncol(geno_complete)
    Hbt <- mean(ifelse(abs(geno_complete[1,]-geno_complete[2,])==1,1,
                       ifelse(geno_complete[1,]==0.5 | geno_complete[2,]==0.5,0.5,0)))
    
    (Hbt-0.5*Hwt-0.5*mean(Hot))/mean(Hot) 
  })
  return(out_rousset_a)
}