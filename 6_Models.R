#### Libs ####

library(dplyr)
library(ape)
library(phytools)
library(nlme)
library(emmeans)
library(circular)
library(mvabund)
library(sandwich)


#### Data input ####

#source("5_F_diversity.R")

tree <- ape::read.tree("Processed data/tree.newick")

model_data <- read.csv(
  "Processed data/imp_data.csv",
  sep = ";",
  dec = ".",
  header = TRUE,
  check.names = FALSE
)

raw_data <- read.csv(
  "Input/Imputation/preimp_data.csv",
  header = TRUE,
  sep = ";",
  dec = ".",
  check.names = FALSE
)

imp_data_nophylo <- read.csv(
  "Processed data/imp_data_nophylo.csv",
  header = TRUE,
  check.names = FALSE
)

#### Settings ####

set.seed(123)

#### Categorical trait groups ####

categorical_groups <- list(
  
  DSyn = c("Autochory", "Ballochory", "Barochory", "Anemochory",
           "Meteorochory", "Endozoochory", "Epizoochory", "Dyszoochory",
           "Myrmecochory", "Hydrochory", "Nautochory"),
  PSyn = c("Anemophily", "Entomophily", "Autogamy"),
  LPhen = c("Evergreen"),
  LShape = c("Acicular", "Cladode", "Claw", "Cordate", "Costapalmate",
             "Cuneate", "Deltoid", "Digitate", "Elliptic", "Falcate",
             "Fan", "Flabellate", "Lanceolate", "Linear", "Lobed",
             "Oblong", "Obovate", "Orbiculate", "Ovate", "Palmate",
             "Phyllode", "Rhomboid", "Scale", "Trifoliate"),
  FlShape = c("Fascicle", "Corymbs", "Actinomorphic", "Glomeroles",
              "Catkins", "Cupule", "Trumpet", "Spadix", "Spike",
              "Compound.dichasium", "Capitulum", "Cymes", "Globose",
              "Raceme", "Filiform", "Panicle", "Umbel"),
  FrShape = c("Cylindrical.Cone", "Pod", "Samara", "Schizocarp",
              "Capsule", "Nut", "Pome", "Berry", "Drupe",
              "Galbulus", "Syncarpus", "Follicle"),
  Btex = c("Ridges", "Smooth", "Lenticels","Peeling.zones", "Peeling.fibrous", "Others"),
  CShape = c("Conical", "Rounded", "Shrub", "Climbing", "Irregular",
             "Fountain", "Columnar", "Weeping", "Prostate", "Layered")
  )

#### Data processing ####

model_data <- model_data |>
  dplyr::mutate(Origin = factor(Origin, levels = c("Native", "Alien")),
    Inv.Status = factor(Inv.Status, levels = c("Native", "Non-Established",
                                               "Naturalised", "Invasive"))
  )

raw_data <- raw_data |>
  dplyr::mutate(Origin = factor(Origin, levels = c("Native", "Alien")),
    Inv.Status = factor(Inv.Status, levels = c("Native", "Non-Established",
                                               "Naturalised", "Invasive"))
  )

model_IS <- model_data |>
  dplyr::filter(Origin == "Alien")

# Datasets for lambda calculation

imp_nophylo_IS <- imp_data_nophylo |>
  dplyr::filter(Origin == "Alien")

#### Pagel's lambda ####

lambda_OR <- list()
lambda_IS <- list()

for (trait in c("SLA", "LNitr", "Dweight", "Height", "Longevity",
                "Tfrost", "Tshade", "Tdrought",
                "minHardiness", "maxHardiness")) {
  
  # Origin (all species)
  
  lambda_data <- imp_data_nophylo |>
    dplyr::transmute(
      Species,
      TraitValue = as.numeric(as.character(.data[[trait]]))
    ) |>
    dplyr::filter(
      is.finite(TraitValue),
      Species %in% tree$tip.label
    )
  
  lambda_tree <- ape::keep.tip(tree, lambda_data$Species)
  
  lambda_data <- lambda_data[
    match(lambda_tree$tip.label, lambda_data$Species),
  ]
  
  values <- lambda_data$TraitValue
  names(values) <- lambda_data$Species
  
  lambda_model <- phytools::phylosig(
    lambda_tree,
    values,
    method = "lambda",
    test = TRUE
  )
  
  lambda_OR[[trait]] <- data.frame(
    Trait = trait,
    N = length(values),
    Lambda = unname(lambda_model$lambda),
    p = unname(lambda_model$P)
  )
  
  # Invasion status
  
  lambda_data_IS <- imp_nophylo_IS |>
    dplyr::transmute(
      Species,
      TraitValue = as.numeric(as.character(.data[[trait]]))
    ) |>
    dplyr::filter(
      is.finite(TraitValue),
      Species %in% tree$tip.label
    )
  
  lambda_tree_IS <- ape::keep.tip(
    tree,
    lambda_data_IS$Species
  )
  
  lambda_data_IS <- lambda_data_IS[
    match(lambda_tree_IS$tip.label, lambda_data_IS$Species),
  ]
  
  values_IS <- lambda_data_IS$TraitValue
  names(values_IS) <- lambda_data_IS$Species
  
  lambda_model_IS <- phytools::phylosig(
    lambda_tree_IS,
    values_IS,
    method = "lambda",
    test = TRUE
  )
  
  lambda_IS[[trait]] <- data.frame(
    Trait = trait,
    N = length(values_IS),
    Lambda = unname(lambda_model_IS$lambda),
    p = unname(lambda_model_IS$P)
  )
}

lambda_OR <- dplyr::bind_rows(lambda_OR)

lambda_IS <- dplyr::bind_rows(lambda_IS)

#### Phylogenetic GLS ####

gls_origin_models <- list()

gls_status_models <- list()

gls_origin_results <- list()

gls_status_contrasts <- list()


for (trait in c("SLA", "LNitr", "Dweight", "Height", "Longevity", "Tfrost", "Tshade", "Tdrought",
                "minHardiness", "maxHardiness")) {
  
  # Origin
  
  model_df <- model_data |>
    dplyr::transmute(Species, Origin, TraitValue = as.numeric(as.character(.data[[trait]]))) |>
    dplyr::filter(is.finite(TraitValue))
  
  model_tree <- ape::keep.tip(tree, model_df$Species)
  
  model_df <- model_df[match(model_tree$tip.label, model_df$Species),]
  
  lambda_value <- lambda_OR |>
    dplyr::filter(Trait == trait) |>
    dplyr::pull(Lambda)
  
  model_origin <- nlme::gls(TraitValue ~ Origin, data = model_df,
                            correlation = ape::corPagel(value = lambda_value, form = ~ Species,
                                                        phy = model_tree, fixed = TRUE),method = "ML")
  
  gls_origin_models[[trait]] <- model_origin
  
  model_summary <- summary(model_origin)$tTable
  
  gls_origin_results[[trait]] <- data.frame(
    Trait = trait,
    Estimate = model_summary[2, "Value"],
    SE = model_summary[2, "Std.Error"],
    t = model_summary[2, "t-value"],
    p = model_summary[2, "p-value"]
    )
  
  # Invasion status
  
  model_df <- model_IS |>
    dplyr::transmute(Species, Inv.Status, TraitValue = as.numeric(as.character(.data[[trait]]))) |>
    dplyr::filter(is.finite(TraitValue))
  
  model_tree <- ape::keep.tip(tree, model_df$Species)
  
  model_df <- model_df[match(model_tree$tip.label, model_df$Species),]
  
  lambda_value <- lambda_IS |>
    dplyr::filter(Trait == trait) |>
    dplyr::pull(Lambda)
  
  model_status <- nlme::gls(TraitValue ~ Inv.Status, data = model_df, 
                            correlation = ape::corPagel(
                              value = lambda_value,form = ~ Species, phy = model_tree, fixed = TRUE),method = "ML")
  
  gls_status_models[[trait]] <- model_status
  
  contrasts <- emmeans::emmeans(model_status, pairwise ~ Inv.Status)$contrasts |>
    summary()
  
  gls_status_contrasts[[trait]] <- data.frame(
    Trait = trait,
    Comparison = contrasts$contrast,
    Estimate = contrasts$estimate,
    SE = contrasts$SE,
    t = contrasts$t.ratio,
    p = contrasts$p.value
    )
}


gls_origin_results <- dplyr::bind_rows(gls_origin_results) |>
  dplyr::mutate(Signif = dplyr::case_when(p < 0.001 ~ "***",
                                          p < 0.01 ~ "**",
                                          p < 0.05 ~ "*",
                                          TRUE ~ ""))

gls_status_contrasts <- dplyr::bind_rows(gls_status_contrasts) |>
  dplyr::mutate(Signif = dplyr::case_when(p < 0.001 ~ "***",
                                          p < 0.01 ~ "**",
                                          p < 0.05 ~ "*",
                                          TRUE ~ ""))

#### Flowering phenology ####

phenology_origin <- list()

phenology_status <- list()

phenology_status_pairs <- list()

for (trait in c("Phen_Time_Start", "Phen_Time_End")) {
  
  # Origin
  
  origin_data <- model_data |>
    dplyr::filter(!is.na(.data[[trait]]))
  
  origin_angles <- circular::circular(origin_data[[trait]] * 2 * pi / 12, units = "radians")
  
  origin_test <- circular::watson.wheeler.test(origin_angles, origin_data$Origin)
  
  phenology_origin[[trait]] <- data.frame(
      Trait = trait,
      Statistic = unname(origin_test$statistic),
      df = unname(origin_test$parameter),
      p = origin_test$p.value
    )
  
  # Invasion status
  
  status_data <- model_IS |>
    dplyr::filter(!is.na(.data[[trait]]))
  
  status_angles <- circular::circular(status_data[[trait]] * 2 * pi / 12,
                                      units = "radians")
  
  status_test <- circular::watson.wheeler.test(status_angles, status_data$Inv.Status)
  
  phenology_status[[trait]] <- data.frame(
      Trait = trait,
      Statistic = unname(status_test$statistic),
      df = unname(status_test$parameter),
      p = status_test$p.value
    )
  
  # Pairwise invasion status
  
  for (pair in list(c("Non-Established","Naturalised"),
                    c("Non-Established", "Invasive"),
                    c("Naturalised", "Invasive"))) {
    
    pair_data <- status_data |>
      dplyr::filter(Inv.Status %in% pair)
    
    pair_angles <- circular::circular(pair_data[[trait]] * 2 * pi / 12,
                                      units = "radians")
    
    pair_test <- circular::watson.wheeler.test(pair_angles, droplevels(pair_data$Inv.Status))
    
    phenology_status_pairs[[paste(trait, pair, collapse = "_")]] <- data.frame(
        Trait = trait,
        Comparison = paste(pair, collapse = " - "),
        Statistic = unname(pair_test$statistic),
        df = unname(pair_test$parameter),
        p = pair_test$p.value
      )
  }
}

phenology_origin <- dplyr::bind_rows(phenology_origin)

phenology_status <- dplyr::bind_rows(phenology_status)

phenology_status_pairs <- dplyr::bind_rows(phenology_status_pairs) |>
  dplyr::group_by(Trait) |>
  dplyr::mutate(p = p.adjust(p, method = "holm")) |>
  dplyr::ungroup() |>
  dplyr::mutate(Signif = dplyr::case_when(p < 0.001 ~ "***",
                                          p < 0.01 ~ "**",
                                          p < 0.05 ~ "*",
                                          TRUE ~ ""))

#### Categorical multivariate models ####

categorical_origin_models <- list()

categorical_status_models <- list()

categorical_origin_anova <- list()

categorical_status_anova <- list()


for (trait_group in names(categorical_groups)) {
  
  traits <- categorical_groups[[trait_group]]
  
  # Origin
  
  cat_origin <- model_data |>
    dplyr::select(Origin, dplyr::all_of(traits))
  
  cat_origin <- cat_origin[complete.cases(cat_origin),]
  
  model_origin <- mvabund::manyglm(mvabund::mvabund(cat_origin[traits]) ~ Origin,
                                   data = cat_origin, family = "binomial")
  
  categorical_origin_models[[trait_group]] <- model_origin
  
  categorical_origin_anova[[trait_group]] <- mvabund::anova.manyglm(model_origin,
                                                                    test = "LR",
                                                                    p.uni = "adjusted")
  
  # Invasion status
  
  cat_status <- model_IS |>
    dplyr::select(Inv.Status, dplyr::all_of(traits))
  
  cat_status <- cat_status[complete.cases(cat_status),]
  
  model_status <- mvabund::manyglm(mvabund::mvabund(cat_status[traits]) ~ Inv.Status,
                                   data = cat_status,family = "binomial")
  
  categorical_status_models[[trait_group]] <- model_status
  
  categorical_status_anova[[trait_group]] <- mvabund::anova.manyglm(model_status,
                                                                    test = "LR",
                                                                    p.uni = "adjusted")
  
  }

#### Categorical individual contrasts ####

categorical_origin_contrasts <- list()

categorical_status_contrasts <- list()


for (trait_group in names(categorical_groups)) {
  
  for (trait in categorical_groups[[trait_group]]) {
    
    # Origin
    
    trait_data <- model_data |>
      dplyr::transmute(Origin, Presence = .data[[trait]])
    
    model_origin <- stats::glm(Presence ~ Origin, data = trait_data,
                               family = stats::binomial())
    
    robust_vcov <- sandwich::vcovHC(model_origin, type = "HC3")
    
    contrast_origin <- emmeans::emmeans(model_origin, pairwise ~ Origin, vcov. = robust_vcov,
                                        adjust = "none")$contrasts |>
      summary()
    
    categorical_origin_contrasts[[paste(trait_group,trait,sep = "_")]] <- data.frame(
        TraitGroup = trait_group,
        Trait = trait,
        Comparison = contrast_origin$contrast,
        Estimate = contrast_origin$estimate,
        SE = contrast_origin$SE,
        z = contrast_origin$z.ratio,
        p = contrast_origin$p.value,
        OR = exp(contrast_origin$estimate),
        OR_low = exp(contrast_origin$estimate - 1.96 * contrast_origin$SE),
        OR_high = exp(contrast_origin$estimate + 1.96 * contrast_origin$SE)
      )
    
    # Invasion status
    
    trait_data <- model_IS |>
      dplyr::transmute(Inv.Status,Presence = .data[[trait]])
    
    model_status <- stats::glm(Presence ~ Inv.Status, data = trait_data, family = stats::binomial())
    
    robust_vcov <- sandwich::vcovHC(model_status, type = "HC3")
    
    contrast_status <- emmeans::emmeans(model_status, pairwise ~ Inv.Status, vcov. = robust_vcov,
                                        adjust = "holm")$contrasts |>
      summary()
    
    categorical_status_contrasts[[paste(trait_group, trait, sep = "_")]] <- data.frame(
      TraitGroup = trait_group,
      Trait = trait,Comparison = contrast_status$contrast,
      Estimate = contrast_status$estimate,
      SE = contrast_status$SE,
      z = contrast_status$z.ratio,
      p = contrast_status$p.value,
      OR = exp(contrast_status$estimate),
      OR_low = exp(contrast_status$estimate - 1.96 * contrast_status$SE),
      OR_high = exp(contrast_status$estimate + 1.96 * contrast_status$SE)
      )
  }
}

categorical_origin_contrasts <- dplyr::bind_rows(categorical_origin_contrasts) |>
  dplyr::mutate(Signif = dplyr::case_when(p < 0.001 ~ "***",
                                          p < 0.01 ~ "**",
                                          p < 0.05 ~ "*",
                                          TRUE ~ "")
                )

categorical_status_contrasts <- dplyr::bind_rows(categorical_status_contrasts) |>
  dplyr::mutate(Signif = dplyr::case_when(p < 0.001 ~ "***",
                                          p < 0.01 ~ "**",
                                          p < 0.05 ~ "*",
                                          TRUE ~ "")
                )

#### Trait amplitude models ####

  # Phen Range

phen_range_origin <- stats::glm(Phen_Range ~ Origin, data = model_data,family = stats::poisson())

phen_range_status <- stats::glm(Phen_Range ~ Inv.Status, data = model_IS,family = stats::poisson())

amplitude_dispersion <- list(
  Phen_Range_Origin = data.frame(
    Trait = "Phen_Range",
    Classification = "Origin",
    Dispersion = sum(residuals(phen_range_origin, type = "pearson")^2) / phen_range_origin$df.residual),
  Phen_Range_Status = data.frame(
    Trait = "Phen_Range",
    Classification = "Inv.Status",
    Dispersion =sum(residuals(phen_range_status, type = "pearson")^2) / phen_range_status$df.residual))

phen_contrasts <- emmeans::emmeans(phen_range_status, pairwise ~ Inv.Status, adjust = "holm")$contrasts |>
  summary()

amplitude_status_contrasts <- list(
  Phen_Range = data.frame(
    Trait = "Phen_Range",
    Comparison = phen_contrasts$contrast,
    Estimate = phen_contrasts$estimate,
    SE = phen_contrasts$SE,
    Statistic = phen_contrasts$z.ratio,
    p = phen_contrasts$p.value
    )
)

  # DSyn & PSyn

amplitude_origin_models <- list()

amplitude_status_models <- list()


for (trait in c("DS_amp", "PS_amp", "LS_amp")) {
  
  poisson_origin <- stats::glm(stats::reformulate("Origin", response = trait),
                               data = model_data,family = stats::poisson())
  
  poisson_status <- stats::glm(stats::reformulate("Inv.Status", response = trait),
                               data = model_IS,family = stats::poisson())
  
  amplitude_dispersion[[paste(trait, "Origin", sep = "_")]] <- data.frame(
    Trait = trait,
    Classification = "Origin",
    Dispersion = sum(residuals(poisson_origin, type = "pearson")^2) / poisson_origin$df.residual
    )
  
  amplitude_dispersion[[paste(trait, "Status", sep = "_")]] <- data.frame(
    Trait = trait,
    Classification = "Inv.Status",
    Dispersion =sum(residuals(poisson_status, type = "pearson")^2) / poisson_status$df.residual
    )
  
  model_origin <- stats::glm(stats::reformulate("Origin", response = trait),
                             data = model_data,family = stats::quasipoisson()
                             )
  
  model_status <- stats::glm(stats::reformulate("Inv.Status", response = trait),
                             data = model_IS,family = stats::quasipoisson()
                             )
  
  amplitude_origin_models[[trait]] <- model_origin
  
  amplitude_status_models[[trait]] <- model_status
  
  contrasts <- emmeans::emmeans(model_status, pairwise ~ Inv.Status, adjust = "holm")$contrasts |>
    summary()
  
  amplitude_status_contrasts[[trait]] <- data.frame(
      Trait = trait,
      Comparison = contrasts$contrast,
      Estimate = contrasts$estimate,
      SE = contrasts$SE,
      Statistic = contrasts$t.ratio,
      p = contrasts$p.value
    )
}


#### Amplitude results ####

amplitude_dispersion <- dplyr::bind_rows(amplitude_dispersion)

amplitude_status_contrasts <- dplyr::bind_rows(amplitude_status_contrasts) |>
  dplyr::mutate(Signif = dplyr::case_when(p < 0.001 ~ "***",
                                          p < 0.01 ~ "**",
                                          p < 0.05 ~ "*",
                                          TRUE ~ "")
  )