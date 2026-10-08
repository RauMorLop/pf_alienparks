#### Libs ####

library(dplyr)
library(ape)
library(circular)
library(gawdis)
library(vegan)
library(TPD)


#### Data input ####

#source("4_F_diversity_num.R")

tree <- ape::read.tree("Processed data/tree.newick")

f_data <- read.csv("Processed data/imp_data.csv", sep = ";", dec = ".",header = TRUE)


#### Settings ####

set.seed(123)

perm <- 999

alpha_tpd <- 0.9

n_divisions_tpd <- 200


#### Data processing ####

f_data <- f_data |>
  dplyr::mutate(Origin = factor(Origin, levels = c("Native", "Alien")),
    Inv.Status = factor(Inv.Status, levels = c("Native", "Non-Established",
                                               "Naturalised", "Invasive"))
  )

#### Correlated traits ####

categorical_traits <- names(f_data)[
  sapply(
    f_data,
    function(x) {
      x <- x[!is.na(x)]
      is.numeric(x) &&
        length(x) > 0 &&
        all(x %in% c(0, 1))
    }
  )
]

cor_mat <- stats::cor(
  f_data[, categorical_traits],
  use = "pairwise.complete.obs"
)

idx <- which(
  abs(cor_mat) >= 0.7 &
    upper.tri(cor_mat),
  arr.ind = TRUE
)

cor_pairs <- data.frame(
  Trait1 = rownames(cor_mat)[idx[, 1]],
  Trait2 = colnames(cor_mat)[idx[, 2]],
  Correlation = cor_mat[idx]
) |>
  dplyr::arrange(
    dplyr::desc(abs(Correlation))
  )

cor_pairs

RP_dispersal <- c("Autochory", "Barochory", "Ballochory", "Anemochory", "Dyszoochory", 
                  "Epizoochory", "Myrmecochory", "Meteorochory", "Hydrochory", "Nautochory")

RP_pollination <- c("Anemophily", "Autogamy")

AES_leaf_phenology <- c("Decidious")

AES_leaf_shape <- c("Acicular", "Cladode", "Claw", "Cordate", "Costapalmate",
                    "Cuneate", "Deltoid", "Digitate", "Elliptic", "Falcate",
                    "Fan", "Flabellate", "Lanceolate", "Linear", "Lobed",
                    "Oblong", "Obovate", "Orbiculate", "Ovate", "Palmate",
                    "Phyllode", "Rhomboid", "Scale", "Trifoliate")

AES_flower_shape <- c("Fascicle", "Corymbs", "Actinomorphic", "Glomeroles", "Catkins",
                      "Cupule", "Trumpet", "Spadix", "Spike", "Compound.dichasium",
                      "Capitulum", "Cymes", "Globose", "Raceme", "Filiform",
                      "Panicle", "Umbel")

AES_fruit_shape <- c("Cylindrical.Cone", "Pod", "Samara", "Schizocarp", "Capsule",
                     "Nut", "Pome", "Berry", "Drupe", "Galbulus",
                     "Syncarpus", "Follicle")

AES_bark_texture <- c("Ridges", "Smooth", "Lenticels", "Peeling.zones", "Peeling.fibrous",
                      "Others")

AES_canopy_shape <- c("Conical", "Rounded", "Shrub", "Climbing", "Irregular",
                      "Columnar", "Weeping", "Prostate","Layered")

#### Data processing ####

 # Circular phenology

f_data <- f_data |>
  dplyr::mutate(
    Phen_Start_rad = 2 * pi * as.numeric(as.character(Phen_Time_Start)) / 12,
    Phen_End_rad = 2 * pi * as.numeric(as.character(Phen_Time_End)) / 12
  )

 # Align data and phylogeny 

common_species <- intersect(tree$tip.label,f_data$Species)

tree <- ape::keep.tip(tree,common_species)

f_data <- f_data[match(tree$tip.label,f_data$Species), ]

 # Reproductive traits

RP_traits <- f_data[, c("Phen_Start_rad", "Phen_End_rad", 
                                  RP_dispersal, RP_pollination),
                              drop = FALSE
                              ]


rownames(RP_traits) <- f_data$Species

RP_traits$Phen_Start_rad <- circular::circular(RP_traits$Phen_Start_rad, units = "radians", modulo = "2pi")

RP_traits$Phen_End_rad <- circular::circular(RP_traits$Phen_End_rad, units = "radians", modulo = "2pi")

RP_groups <- c(
  rep("Phenology", 2),
  rep("Dispersal Syndrome", length(RP_dispersal)),
  rep("Pollination Syndrome", length(RP_pollination))
)

names(RP_groups) <- colnames(RP_traits)

RP_fuzzy <- c(
  "Dispersal Syndrome",
  "Pollination Syndrome"
)

  # Aesthetic traits

AES_traits <- f_data[,c(AES_leaf_phenology, AES_leaf_shape, AES_flower_shape,
                                  AES_fruit_shape, AES_bark_texture, AES_canopy_shape),
                               drop = FALSE]

rownames(AES_traits) <- f_data$Species

AES_groups <- c(
  rep("Leaf Phenology", length(AES_leaf_phenology)),
  rep("Leaf Shape", length(AES_leaf_shape)),
  rep("Flower Shape", length(AES_flower_shape)),
  rep("Fruit Shape", length(AES_fruit_shape)),
  rep("Bark Texture", length(AES_bark_texture)),
  rep("Canopy Shape",length(AES_canopy_shape))
)

names(AES_groups) <- colnames(AES_traits)

AES_fuzzy <- c(
  "Leaf Shape",
  "Flower Shape",
  "Fruit Shape",
  "Bark Texture",
  "Canopy Shape"
)

 # Trait groups and fuzzy traits

trait_data <- list(RP = RP_traits, AES = AES_traits)

trait_groups <- list(RP = RP_groups, AES = AES_groups)

trait_fuzzy <- list(RP = RP_fuzzy, AES = AES_fuzzy)

#### IMPORT REND_df() from the "Functions" folder ####

#### distance matrices, TPD and functional diversity ####

distance_matrices <- list()

pcoa_objects <- list()

tpd_objects <- list()

ordination_results <- list()

vector_results <- list()

variance_results <- list()

mantel_results <- list()

FD_results <- list()


for (trait_group in names(trait_data)) {
  
  traits <- trait_data[[trait_group]]
  
  groups <- trait_groups[[trait_group]]
  
  fuzzy <- trait_fuzzy[[trait_group]]
  
  # Functional distance matrix 
  
  FDist <- gawdis::gawdis(traits, w.type = "optimized", groups = groups, 
                          groups.weight = FALSE, fuzzy = fuzzy)
  
  distance_matrices[[trait_group]] <- FDist
  
  # Mantel test
  
  species <- attr(FDist, "Labels")
  
  PDist <- ape::cophenetic.phylo(tree)
  
  PDist <- PDist[species, species, drop = FALSE]
  
  mantel_model <- vegan::mantel(FDist, stats::as.dist(PDist), method = "pearson",
                                permutations = perm)
  
  mantel_results[[trait_group]] <- data.frame(
      TraitGroup = trait_group,
      Mantel_r = as.numeric(mantel_model$statistic),
      p = mantel_model$signif
    )
  
  
  # PCoA 
  
  pcoa <- vegan::wcmdscale(FDist, eig = TRUE, add = FALSE)
  
  pcoa_objects[[trait_group]] <- pcoa
  
  eig <- pcoa$eig[pcoa$eig > 0]
  
  variance_results[[trait_group]] <- data.frame(
      TraitGroup = trait_group,
      Axis = seq_along(eig),
      Variance = eig * 100 / sum(eig)
    )
  
  coords <- as.data.frame(pcoa$points[,1:2,drop = FALSE])
  
  colnames(coords) <- c("PCo1", "PCo2")
  
  ordination_results[[trait_group]] <- data.frame(
    TraitGroup = trait_group,
    Species = rownames(coords),
    Family = f_data$Family,
    Origin = f_data$Origin,
    Inv.Status = f_data$Inv.Status,
    PCo1 = coords$PCo1,
    PCo2 = coords$PCo2
    )
  
  # Fitted vectors
  
  traits_env <- as.data.frame(lapply(traits,function(x) {
        
        if (inherits(x,"circular")) {
          as.numeric(x)
        } else {
          as.numeric(as.character(x))
        }
      }
    ))
  
  rownames(traits_env) <- rownames(traits)
  
  fit <- vegan::envfit(as.matrix(coords), traits_env, permutations = 0)
  
  arrows <- as.data.frame(fit$vectors$arrows)
  
  correlations <- stats::cor(traits_env, as.matrix(coords), use = "pairwise.complete.obs")
  
  vector_results[[trait_group]] <- data.frame(
    TraitGroup = trait_group,
    Trait = rownames(arrows),
    PCo1 = arrows[, 1],
    PCo2 = arrows[, 2],
    R2 = fit$vectors$r[rownames(arrows)],
    Cor_PCo1 = correlations[rownames(arrows), 1],
    Cor_PCo2 = correlations[rownames(arrows),2]
    )
  
  # TPD for Origin grops
  
  coords_matrix <- as.matrix(coords)
  
  rownames(coords_matrix) <- rownames(coords)
  
  TPD_origin <- TPD::TPDs(species = as.character(f_data$Origin),
                          traits = coords_matrix,
                          alpha = alpha_tpd,
                          n_divisions = n_divisions_tpd)
  
  # Spatial grid
  
  grid_origin <- as.data.frame(TPD_origin$data$evaluation_grid)
  
  trait_ranges <- list(
    range(grid_origin[[1]],na.rm = TRUE),
    range(grid_origin[[2]],na.rm = TRUE))
  
  n_divisions <- length(unique(grid_origin[[1]]))
  
  # TPD for Invasion status groups
  
  alien_idx <- which(f_data$Origin =="Alien" & !is.na(f_data$Inv.Status))
  
  TPD_status <- TPD::TPDs(
    species = as.character(f_data$Inv.Status[alien_idx]),
    traits = coords_matrix[alien_idx, ,drop = FALSE],
    alpha = alpha_tpd,
    trait_ranges = trait_ranges,
    n_divisions = n_divisions
    )
  
  tpd_objects[[trait_group]] <- list(Origin = TPD_origin, Inv.Status = TPD_status)
  
  # Functional diversity indices and null models
  
  for (classification in c("Origin","Inv.Status")) {
    
    if (classification == "Origin") {
      
      TPD_obs <- TPD_origin
      
      labels <- as.character(f_data$Origin)
      
      traits_tpd <- coords_matrix
      
    } else {
      
      TPD_obs <- TPD_status
      
      labels <- as.character(f_data$Inv.Status[alien_idx])
      
      traits_tpd <- coords_matrix[alien_idx, , drop = FALSE]
    }
    
    # Obs. Functional indices
    
    obs <- REND_df(TPD::REND(TPDs = TPD_obs))
    
    # Grid
    
    grid <- as.data.frame(TPD_obs$data$evaluation_grid)
    
    null_ranges <- list(range(grid[[1]], na.rm = TRUE), range(grid[[2]], na.rm = TRUE))
    
    null_divisions <- length(unique(grid[[1]]))
    
    # Null models
    
    null_list <- vector("list",perm)
    
    for (perm in seq_len(perm)) {
      
      TPD_null <- TPD::TPDs(species = sample(labels),
                            traits = traits_tpd,
                            alpha = alpha_tpd,
                            trait_ranges = null_ranges,
                            n_divisions = null_divisions)
      
      null_list[[perm]] <- REND_df(TPD::REND(TPDs = TPD_null)) |>
        dplyr::mutate(Sim = perm)
    }
    
    null_df <- dplyr::bind_rows(null_list)
    
    
    # Null distributions
    
    null_summary <- null_df |>
      dplyr::group_by(Group, Metric) |>
      dplyr::summarise(Mean_null = mean(Obs, na.rm = TRUE),
                       SD_null = sd(Obs, na.rm = TRUE),
                       Null_values = list(Obs[!is.na(Obs)]),
                       .groups = "drop")
    
    # SES and p-values 
    
    FD_results[[paste(trait_group, classification, sep = "_")]] <- obs |>
      dplyr::left_join(null_summary, by = c("Group", "Metric")) |>
      dplyr::mutate(SES = ifelse(
        is.na(SD_null) | SD_null == 0, NA_real_,
        (Obs - Mean_null) / SD_null),
        P = mapply(
          
          function(null,mu,obs) {
            (sum(abs(null - mu) >= abs(obs - mu)) + 1) / (length(null) + 1)
            },
          Null_values,
          Mean_null,
          Obs
          ),
        Signif = dplyr::case_when(P < 0.001 ~ "***", 
                                  P < 0.01 ~ "**", 
                                  P < 0.05 ~ "*", 
                                  TRUE ~ ""),
        TraitGroup = trait_group,
        Classification = classification,
        .before = 1
      ) |>
      dplyr::select(-Null_values)
  }
}


#### Results ####

ordination_results <- dplyr::bind_rows(ordination_results)

vector_results <- dplyr::bind_rows(vector_results)

variance_results <- dplyr::bind_rows(variance_results)

mantel_results <- dplyr::bind_rows(mantel_results)

FD_results <- dplyr::bind_rows(FD_results)