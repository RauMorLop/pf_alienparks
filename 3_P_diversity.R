#### Libs ####

library(ape)
library(picante)
library(dplyr)


#### Data input ####

#source("2_Imputation.R")

tree <- ape::read.tree("Processed data/tree.newick")

P_Dist <- ape::cophenetic.phylo(tree)

phyl_data <- read.csv("Processed data/imp_data.csv", header = TRUE, sep = ";", dec = ".")

runs <- 999

set.seed(123)


#### Data processing ####

all_tips <- tree$tip.label

groups <- list(
  Native = phyl_data |>
    filter(Origin == "Native") |>
    pull(Species),
  
  Alien = phyl_data |>
    filter(Origin == "Alien") |>
    pull(Species),
  
  `Non-Established` = phyl_data |>
    filter(Inv.Status == "Non-Established") |>
    pull(Species),
  
  Naturalised = phyl_data |>
    filter(Inv.Status == "Naturalised") |>
    pull(Species),
  
  Invasive = phyl_data |>
    filter(Inv.Status == "Invasive") |>
    pull(Species)
)

#### Phylogenetic diversity ####

phyl_results <- list()

for (group in names(groups)) {
  
  spp <- intersect(
    groups[[group]],
    all_tips
  )
  
  comm <- matrix(
    0,
    nrow = 1,
    ncol = length(all_tips),
    dimnames = list(group, all_tips)
  )
  
  comm[1, spp] <- 1
  
  
  # Duplicate row for picante SES functions
  
  comm2 <- rbind(comm, comm)
  
  rownames(comm2) <- c(
    group,
    "Dummy"
  )
  
  
  #### Faith's PD ####
  
  PD <- picante::ses.pd(
    comm2,
    tree,
    null.model = "taxa.labels",
    runs = runs
  )[1, ]
  
  
  #### PR and PV ####
  
  PR_obs <- picante::psr(
    as.data.frame(comm),
    tree,
    compute.var = FALSE
  )$PSR[1]
  
  PV_obs <- picante::psv(
    as.data.frame(comm),
    tree,
    compute.var = FALSE
  )$PSV[1]
  
  
  # Null distributions for PR and PV
  
  null_values <- replicate(
    runs,
    {
      
      null_tree <- tree
      
      null_tree$tip.label <- sample(
        null_tree$tip.label
      )
      
      c(
        PR = picante::psr(
          as.data.frame(comm),
          null_tree,
          compute.var = FALSE
        )$PSR[1],
        
        PV = picante::psv(
          as.data.frame(comm),
          null_tree,
          compute.var = FALSE
        )$PSV[1]
      )
    }
  )
  
  
  # PR
  
  PR_null <- null_values["PR", ]
  
  PR_mean <- mean(PR_null)
  PR_sd <- sd(PR_null)
  
  SES_PR <- (
    PR_obs - PR_mean
  ) / PR_sd
  
  PR_p <- (
    sum(
      abs(PR_null - PR_mean) >=
        abs(PR_obs - PR_mean)
    ) + 1
  ) / (runs + 1)
  
  
  # PV
  
  PV_null <- null_values["PV", ]
  
  PV_mean <- mean(PV_null)
  PV_sd <- sd(PV_null)
  
  SES_PV <- (
    PV_obs - PV_mean
  ) / PV_sd
  
  PV_p <- (
    sum(
      abs(PV_null - PV_mean) >=
        abs(PV_obs - PV_mean)
    ) + 1
  ) / (runs + 1)
  
  
  #### NRI ####
  
  MPD <- picante::ses.mpd(
    comm2,
    P_Dist,
    null.model = "taxa.labels",
    abundance.weighted = FALSE,
    runs = runs
  )[1, ]
  
  
  #### NTI ####
  
  MNTD <- picante::ses.mntd(
    comm2,
    P_Dist,
    null.model = "taxa.labels",
    abundance.weighted = FALSE,
    runs = runs
  )[1, ]
  
  
  #### Results ####
  
  phyl_results[[group]] <- data.frame(
    Group = group,
    S = length(spp),
    
    PD = PD$pd.obs,
    SES_PD = PD$pd.obs.z,
    PD_p = PD$pd.obs.p,
    
    PR = PR_obs,
    SES_PR = SES_PR,
    PR_p = PR_p,
    
    PV = PV_obs,
    SES_PV = SES_PV,
    PV_p = PV_p,
    
    NRI = -MPD$mpd.obs.z,
    NRI_p = MPD$mpd.obs.p,
    
    NTI = -MNTD$mntd.obs.z,
    NTI_p = MNTD$mntd.obs.p
  )
}


#### Results ####

phyl_results <- dplyr::bind_rows(
  phyl_results
)

phyl_results
