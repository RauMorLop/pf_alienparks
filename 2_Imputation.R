#### Libs ####

library(missForest)
library(ape)
library(dplyr)
library(WorldFlora)


#### Data input ####

source("1_Phylogenetic_tree.R")

preimp_data <- read.csv(
  "Input/Imputation/preimp_data.csv",
  header = TRUE,
  dec = ".",
  sep = ";"
)

metadata <- preimp_data |>
  dplyr::select(Species, Family, Origin, Inv.Status)

log_traits <- c("SLA", "LNitr", "Dweight", "Height", "Longevity")

phenology_traits <- c("Phen_Time_Start", "Phen_Time_End")

hardiness_traits <- c("minHardiness", "maxHardiness")

non_dummy_traits <- c(log_traits, phenology_traits, "Phen_Range", "DS_amp", "PS_amp", "Tfrost",
                      "Tshade", "Tdrought", hardiness_traits, "HardRange", "BT_amp", "Flsh_amp", "LS_amp")

dummy_traits <- setdiff(names(preimp_data), c("Species", "Family", "Origin", "Inv.Status", non_dummy_traits))

#### Data transformation ####

  # Histograms to check whether numerical data needs log-transformation

hist(preimp_data$Dweight)
hist(preimp_data$Height)
hist(preimp_data$Longevity)
hist(preimp_data$SLA)
hist(preimp_data$LNitr)

transf_preimp_data <- preimp_data |>
  dplyr::select(-Family, -Origin, -Inv.Status) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(phenology_traits), ~ factor(.x, levels = 1:12, ordered = FALSE)),
                dplyr::across(dplyr::all_of(hardiness_traits), ~ factor(.x, levels = 1:13, ordered = TRUE)),
                dplyr::across(dplyr::all_of(dummy_traits), ~ factor(.x, levels = c(0, 1), ordered = FALSE)),
                Dweight = as.numeric(Dweight)) |>
  dplyr::mutate(dplyr::across(dplyr::all_of(log_traits), ~ log(.x)))

#### Species names match ####

#### Phylogenetic PCoA ####

P_Dist <- as.matrix(P_Dist)

pPCoA <- ape::pcoa(P_Dist)

pPCoA$vectors

  # PCoA visualization

plot(
  pPCoA$vectors[, 1],
  pPCoA$vectors[, 2],
  xlab = "PCo 1",
  ylab = "PCo 2",
  main = "PCoA",
  pch = 19
)

pcoa_eigval <- pPCoA$values$Relative_eig

barplot(
  pcoa_eigval,
  main = "Variance explained"
)

sum(pcoa_eigval[1:4])

  # Fitted vectors

pcoa_vectors <- as.data.frame(pPCoA$vectors[, 1:4])

colnames(pcoa_vectors) <- paste0("PCoA", seq_len(ncol(pcoa_vectors)))

pcoa_vectors$Species <- rownames(pcoa_vectors)

preimp_data_phylo <- transf_preimp_data |>
  dplyr::left_join(pcoa_vectors, by = "Species")

#### Data imputation ####

imputation_data <- preimp_data_phylo |>
  dplyr::select(-Species,-Family,-Origin,-Inv.Status)

set.seed(123)

imputation <- missForest::missForest(
  imputation_data,
  ntree = 300,
  variablewise = TRUE
)

str(imputation)

  # Imputation error

imp_error <- imputation$OOBerror

imp_error

  # Imputed data

imp_data <- imputation$ximp

imp_data <- imputation$ximp |>
  dplyr::select(
    -dplyr::starts_with("PCoA")
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        c(
          phenology_traits,
          hardiness_traits,
          dummy_traits
        )
      ),
      ~ as.numeric(
        as.character(.x)
      )
    ),
    dplyr::across(
      dplyr::all_of(log_traits),
      exp
    )
  )

imp_data <- cbind(
  metadata,
  imp_data
)

#### Data imputation without phylogeny for phyl. signal ####

imp_nophylo_data <- trans_preimp_data |>
  dplyr::select(-Species)

set.seed(123)

imputation_nophylo <- missForest::missForest(
  transf_preimp_data,
  ntree = 300,
  variablewise = TRUE
)

  # Imputation error

imp_error_nophylo <- imputation_nophylo$OOBerror

imp_error_nophylo

# The OOBerror is high, but we found similar phylogenetic signal values
# and the same significant differences in the models when we imputed with
# and without phylogeny vectors for every trait but SLA.
# The OOBerror for SLA imputation without phylogenetic vectors is low, so
# we will keep this imputation for phylogenetic signal calculation.

  # Imputed data 

imp_data_nophylo <- imputation_nophylo$ximp |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(
        c(
          phenology_traits,
          hardiness_traits,
          dummy_traits
        )
      ),
      ~ as.numeric(
        as.character(.x)
      )
    ),
    dplyr::across(
      dplyr::all_of(log_traits),
      exp
    )
  )

imp_data_nophylo <- cbind(
  metadata,
  imp_data_nophylo
)

#### Save processed data ####

write.csv(
  imp_data,
  "Processed data/imp_data.csv",
  row.names = FALSE
)

write.csv(
  imp_data_nophylo,
  "Processed data/imp_data_nophylo.csv",
  row.names = FALSE
)
