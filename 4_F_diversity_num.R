#### Libs ####

library(dplyr)
library(ape)
library(igraph)
library(phytools)
library(TPD)


#### Data input ####

#source("3_P_diversity.R")

tree <- ape::read.tree("Processed data/tree.newick")

f_data <- read.csv("Processed data/imp_data.csv", sep =";", dec = ".", header = TRUE)

#### Settings ####

set.seed(123)

perm <- 999
alpha_tpd <- 0.9
cor_cutoff <- 0.7


#### Data processing ####

f_data <- f_data |>
  dplyr::mutate(
    Origin = factor(Origin, levels = c("Native", "Alien")),
    Inv.Status = factor(Inv.Status, levels = c("Native", "Non-Established",
                                               "Naturalised", "Invasive"))
  )

#### Correlations ####

cor_traits <- list(
  AUR = numeric_traits$AUR,
  TAS = c("Tfrost", "Tshade", "Tdrought", "minHardiness", "maxHardiness")
)

cor_results <- list()

for (group in names(cor_traits)) {
  
  cor_mat <- stats::cor(
    f_data[,cor_traits[[group]],drop = FALSE],
    use = "pairwise.complete.obs"
  )
  
  idx <- which(
    abs(cor_mat) > cor_cutoff &
      upper.tri(cor_mat),
    arr.ind = TRUE
  )
  
  if (nrow(idx) > 0) {
    
    cor_results[[group]] <- data.frame(
      TraitGroup = group,
      Trait1 = rownames(cor_mat)[idx[, 1]],
      Trait2 = colnames(cor_mat)[idx[, 2]],
      Correlation = cor_mat[idx]
      )
  }
}

cor_results <- dplyr::bind_rows(
  cor_results
)

numeric_traits <- list(
  AUR = c("SLA", "LNitr", "Dweight", "Height", "Longevity"),
  TAS = c("Tfrost", "Tshade", "Tdrought", "maxHardiness")
)

#### Numerical trait table ####

trait_names <- unique(
  unlist(
    numeric_traits,
    use.names = FALSE
  )
)

num_data <- f_data |>
  dplyr::select(
    Species,
    Family,
    Origin,
    Inv.Status,
    dplyr::all_of(trait_names)
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        numeric_traits$AUR
      ),
      log
    )
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        trait_names
      ),
      ~ as.numeric(
        scale(.x)
      )
    )
  )

#### Align species names ####

common_sp <- intersect(tree$tip.label,num_data$Species)

tree <- ape::keep.tip(tree,common_sp)

num_data <- num_data[
  match(tree$tip.label,num_data$Species),
]

source("Functions/REND_df.R")

#### p-PCA, TPD and functional diversity ####

pca_objects <- list()

tpd_objects <- list()

overlap_objects <- list()


ordination_results <- list()

loading_results <- list()

variance_results <- list()

FD_results <- list()


for (trait_group in names(numeric_traits)) {

  traits <-
    numeric_traits[[trait_group]]
  
  trait_matrix <- num_data |>
    dplyr::select(dplyr::all_of(traits)) |>
    as.data.frame()
  
  rownames(trait_matrix) <- num_data$Species
  
  # Phylogenetic PCA
  
  pca <- phytools::phyl.pca(tree, trait_matrix, method = "lambda")
  
  pca_objects[[trait_group]] <-pca
  
  # scores
  
  coords <- as.data.frame(pca$S[, 1:2, drop = FALSE])
  
  colnames(coords) <- c("PC1","PC2")
  
  scores <- data.frame(  
    TraitGroup = trait_group,
    Species = rownames(coords),
    Family = num_data$Family,
    Origin = num_data$Origin,
    Inv.Status = num_data$Inv.Status,
    PC1 = coords$PC1,
    PC2 = coords$PC2
  )
  
  ordination_results[[trait_group]] <- scores
  
  
  #### p-PCA loadings ####
  
  loadings <- as.data.frame(pca$L[, 1:2, drop = FALSE])
  
  colnames(loadings) <- c("PC1", "PC2")
  
  loading_results[[trait_group]] <- 
    data.frame(
      TraitGroup = trait_group,
      Trait = rownames(loadings),
      PC1 = loadings$PC1,
      PC2 = loadings$PC2
    )
  
  # Explained variance
  
  variance_results[[trait_group]] <- data.frame(
      TraitGroup = trait_group,
      Axis = seq_along(pca$Eval),
      Variance = pca$Eval *100 / sum(pca$Eval)
  )
  
  # TPD for Origin groups 
  
  coord_matrix <- as.matrix(scores[, c("PC1", "PC2")])
  
  rownames(coord_matrix) <- scores$Species
  
  TPD_origin <- TPD::TPDs(
    species = as.character(scores$Origin),
    traits = coord_matrix,
    alpha = alpha_tpd
  )
  
  # Spatial grid 
  
  grid <- as.data.frame(TPD_origin$data$evaluation_grid)
  
  trait_ranges <- list(
    range(grid[[1]],na.rm = TRUE),
    range(grid[[2]],na.rm = TRUE)
  )
  
  n_divisions <- c(
    length(unique(grid[[1]])),
    length(unique(grid[[2]]))
  )

  # TPD for Invasion status groups
  
  alien_idx <- which(scores$Origin == "Alien" & !is.na(scores$Inv.Status))
  
  TPD_status <- TPD::TPDs(
    species = as.character(scores$Inv.Status[alien_idx]),
    traits = coord_matrix[alien_idx, , drop = FALSE],
    alpha = alpha_tpd,
    trait_ranges = trait_ranges,  
    n_divisions = n_divisions
  )
  
  tpd_objects[[trait_group]] <- list(
    Origin = TPD_origin,
    Inv.Status = TPD_status
    )
  

  # Functional diversity & SES 
  
  for (classification in c("Origin", "Inv.Status")) {
    
    if (classification == "Origin") {
      
      TPD_obs <- TPD_origin
      
      labels <- as.character(scores$Origin)
      
      traits_tpd <- coord_matrix
      
    } else {
      
      TPD_obs <- TPD_status
      
      labels <- as.character(scores$Inv.Status[alien_idx])
      
      traits_tpd <- coord_matrix[alien_idx, ,drop = FALSE]
    }
    
    # Obs. functional indices
    
    obs <- REND_df(TPD::REND(TPDs = TPD_obs))
    
    # Grid 
    
    grid <- as.data.frame(TPD_obs$data$evaluation_grid)
    
    null_ranges <- list(
      range(grid[[1]], na.rm = TRUE),
      range(grid[[2]], na.rm = TRUE)
    )
    
    null_divisions <- c(
      length(unique(grid[[1]])),
      length(unique(grid[[2]]))
    )
    
    # Null models

    null_list <- vector("list",perm)
    
    for (perm in seq_len(perm)) {
      
      TPD_null <- TPD::TPDs(
        species = sample(labels),  
        traits = traits_tpd,  
        alpha = alpha_tpd,  
        trait_ranges = null_ranges,
        n_divisions = null_divisions
      )
      
      null_list[[perm]] <- REND_df(
          TPD::REND(TPDs = TPD_null)) |>
        dplyr::mutate(Sim = perm)
    }
    
    null_df <- dplyr::bind_rows(null_list)
    
    # Null distributions 
    
    null_summary <- null_df |>
      dplyr::group_by(Group,Metric) |>
      dplyr::summarise(
        Mean_null = mean(Obs, na.rm = TRUE),
        SD_null = sd(Obs, na.rm = TRUE),
        N_null = sum(!is.na(Obs)),
        Null_values = list(Obs[!is.na(Obs)]),
        .groups = "drop"
      )
    
    
    #### SES
    
    FD_results[[paste(trait_group,classification,sep = "_")]] <- obs |>
      dplyr::left_join(null_summary,by = c("Group","Metric")) |>
      dplyr::mutate(SES = ifelse(
        is.na(SD_null) | SD_null == 0, NA_real_,
        (Obs - Mean_null) / SD_null),
        p = mapply(
          
          function(null,mu,obs) {
            (sum(abs(null - mu) >= abs(obs - mu)) +1) / (length(null) + 1)
            },
          Null_values,
          Mean_null,
          Obs
          ),
        Signif = dplyr::case_when(p < 0.001 ~ "***", 
                                  p < 0.01 ~ "**", 
                                  p < 0.05 ~ "*", 
                                  TRUE ~ ""),
        TraitGroup = trait_group,
        Classification = classification,
        .before = 1
      ) |>
      dplyr::select(
        -Null_values)
  }
}

#### Results ####

ordination_results <- dplyr::bind_rows(ordination_results)

loading_results <- dplyr::bind_rows(loading_results)

variance_results <- dplyr::bind_rows(variance_results)

FD_results <- dplyr::bind_rows(FD_results)
