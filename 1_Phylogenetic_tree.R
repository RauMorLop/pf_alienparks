
#### Libs ####

#devtools::install_github("jinyizju/V.PhyloMaker2")

library(devtools)
library(V.PhyloMaker2)
library(ape)
library(dplyr)

#### Data input ####

sp_tree <- read.csv("Input/sp_tree.csv", header = T, sep = ";")

#### Species name standardisation ####

sp_tree$Species <- gsub(" ", "_", sp_tree$Species)

sp_tree$genus <- sub("_.*", "", sp_tree$Species) 

#### Phylogenetic tree ####

sp_tree <- sp_tree |>
  dplyr::transmute(
    species = Species,
    genus = genus,
    family = Family
  )

tree <- phylo.maker(
  sp.list = sp_tree,
  tree = GBOTB.extended.TPL,
  scenarios = "S3")$scenario.3

ape::write.tree(
  tree,
  file = "Processed data/tree.newick"
)

plot(tree)

#### Phylogenetic distances ####

P_Dist <- cophenetic.phylo(tree)
