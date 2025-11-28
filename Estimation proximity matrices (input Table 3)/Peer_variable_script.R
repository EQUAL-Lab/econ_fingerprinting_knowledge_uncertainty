# ------------------------------------------------------------------------
# Author: Luis Sanchez (ls2252@cornell.edu)
# Manuscript: Econ fingerprinting
# Description: Script which codes the variable of proximity
# ------------------------------------------------------------------------

# Load packages
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Load directories (when replicating, use your local directory)
path_base  <- "C:/Users/lalfo/OneDrive/Escritorio/Cornell/2025/EconPaperFinger/Nuevo Peer"
path_data  <- file.path(path_base, "fingerprint_updated_data.csv")
path_north <- file.path(path_base, "Matriz_Adjacencia_Norte_Nueva.csv")
path_south <- file.path(path_base, "Matriz_Adjacencia_Sur_Nueva.csv")

df <- read_csv(path_data, show_col_types = FALSE) %>%
  mutate(IDPLOT = as.character(IDPLOT))

# Load the adjacency matrices of the distances between neighbors in KM (these adjacency matrices were calculated in the QGIS program)

read_distance_matrix <- function(path_csv) {
  raw <- read_csv(path_csv, show_col_types = FALSE)
  id_col <- names(raw)[1]
  ids_row <- as.character(raw[[id_col]])
  mat_vals <- as.matrix(raw[,-1, drop = FALSE])
  colnames(mat_vals) <- make.names(colnames(raw)[-1], unique = TRUE)
  
  # Match original name columns
  original_colnames <- colnames(raw)[-1]
  names(original_colnames) <- colnames(mat_vals)
  
  colnames(mat_vals) <- original_colnames
  
  # Double check that columns name correspondend
  if(!all(original_colnames %in% ids_row)){
    warning("Not all column names appear in the first column (row IDs). The intersection will continue.")
  }
  
  common_ids <- intersect(ids_row, original_colnames)
  if(length(common_ids) == 0) stop("There are no common IDs between rows and columns of the matrix.")
  
  row_idx <- match(common_ids, ids_row)
  col_idx <- match(common_ids, original_colnames)
  
  M <- mat_vals[row_idx, col_idx, drop = FALSE]
  dimnames(M) <- list(common_ids, common_ids)
  
  storage.mode(M) <- "numeric"
  
  return(M)
}

# Compute peer variables
compute_peer_vars <- function(M, df, region_name, radii_km = c(1, 5), match_var = "match2") {
  ids_M <- rownames(M)
  
  # Exclude distance from the respondent him/herself 
  diag(M) <- Inf
  
  # Create a "match" vector following ids_M
  # (IDs which are not in the dataset df will be returned as NA)
  match_vec <- df %>%
    select(IDPLOT, all_of(match_var)) %>%
    mutate(IDPLOT = as.character(IDPLOT)) %>%
    right_join(tibble(IDPLOT = ids_M), by = "IDPLOT") %>%
    pull(all_of(match_var))
  
  # We exclude plots which belong to the same farm
  # Vector with farm ID (IDFARM)
  farm_vec <- df %>%
    select(IDPLOT, IDFARM) %>%
    mutate(IDPLOT = as.character(IDPLOT),
           IDFARM  = as.character(IDFARM)) %>%
    right_join(tibble(IDPLOT = ids_M), by = "IDPLOT") %>%
    pull(IDFARM)
  
  # Boolean matrix: TRUE if i y j are from the same farm (NA otherwise)
  same_farm <- outer(farm_vec, farm_vec, function(a, b) !is.na(a) & !is.na(b) & a == b)
  
  M[same_farm] <- Inf
  
  # For each radio, calculate count and share of peer farmers who also score a match
  out_list <- lapply(radii_km, function(rk){
    within_r <- (M <= rk) # logic matix of neighbours within the radius
    
    # Erase the farm of the same respondents
    within_r <- within_r & !same_farm
    
    
    # Count of neighbors in the radius (excluding the diagonal and the main farm of the respondent)
    neigh_count <- rowSums(within_r, na.rm = TRUE)
    
    # Count of neighbors with match = 1
    match_mat <- matrix(match_vec, nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
    count_match <- rowSums(within_r & (match_mat == 1), na.rm = TRUE)
    
    # Share (NA if neigh_count==0)
    share_match <- ifelse(neigh_count > 0, count_match / neigh_count, NA_real_)
    
    tibble(
      IDPLOT = ids_M,
      !!paste0("count_neigh_", rk, "km_", region_name) := as.numeric(neigh_count),
      !!paste0("count_match_", rk, "km_", region_name) := as.numeric(count_match),
      !!paste0("ratio_match_", rk, "km_", region_name) := as.numeric(share_match)
    )
  })
  
  # Merge with ID plot in the main dataset
  region_df <- reduce(out_list, left_join, by = "IDPLOT")
  
  region_df
}

# Read the two matrices 
M_north <- read_distance_matrix(path_north)
M_south <- read_distance_matrix(path_south)

# Calculate by region
peer_north <- compute_peer_vars(M_north, df, region_name = "North", radii_km = c(1,5), match_var = "match2")
peer_south <- compute_peer_vars(M_south, df, region_name = "South", radii_km = c(1,5), match_var = "match2")

# Join the main datasets 
df_out <- df %>%
  left_join(peer_north, by = "IDPLOT") %>%
  left_join(peer_south, by = "IDPLOT")

write_csv(df_out, file.path(path_base, "Data_Peer_with_peerVars.csv"))

print(df_out %>% select(IDPLOT, starts_with("count_match_"), starts_with("share_match_")) %>% head())
