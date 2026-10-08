#### Libs ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#### Data input ####

source("6_Models.R")

#### Settings ####

origin_cols <- c(
  Native = "#2C7FB8",
  Alien = "#D7301F"
)

status_cols <- c(
  `Non-Established` = "#006400",
  Naturalised = "#6A3D9A",
  Invasive = "#FDBF11"
)

heat_cols <- c(
  "#ffff00",
  "#fdae61",
  "#d7191c"
)


#### Functions for TPD heatmaps ####

source("Functions/TPD_to_plot.R")
source("Functions/TPD_thresholds.R")


tpd_panel <- function(tpd, focal, colours, arrows,
                      axes, xlim, ylim, arrow_scale,
                      contour_levels, contour_width) {
  
  df <- tpd_to_plot(tpd)
  
  
  # Density thresholds
  
  thresholds <- df |>
    dplyr::group_by(Group) |>
    dplyr::summarise(
      t95 = tpd_thresholds(Density)["t95"],
      t50 = tpd_thresholds(Density)["t50"],
      t20 = tpd_thresholds(Density)["t20"],
      t05 = tpd_thresholds(Density)["t05"],
      zmax = max(Density, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  # Focal group
  
  focal_df <- df |>
    dplyr::filter(Group == focal) |>
    dplyr::left_join(
      thresholds,
      by = "Group"
    ) |>
    dplyr::mutate(
      z_rel = pmax(
        (Density - t95) /
          pmax(zmax - t95, 1e-12),
        0
      ),
      z_rel = ifelse(
        Density >= t95,
        z_rel,
        NA_real_
      )
    )
  
  
  # Non-focal groups
  
  other_groups <- setdiff(
    unique(df$Group),
    focal
  )
  
  p <- ggplot()
  
  
  # Contours of non-focal groups
  
  for (group in other_groups) {
    
    contour_data <- df[
      df$Group == group,
    ]
    
    group_thresholds <- thresholds[
      thresholds$Group == group,
    ]
    
    for (level in contour_levels) {
      
      p <- p +
        geom_contour(
          data = contour_data,
          aes(
            Axis1,
            Axis2,
            z = Density,
            colour = Group
          ),
          breaks = group_thresholds[[level]],
          linewidth = contour_width,
          show.legend = FALSE
        )
    }
  }
  
  
  # Heatmap of focal group
  
  p <- p +
    geom_raster(
      data = focal_df,
      aes(
        Axis1,
        Axis2,
        fill = z_rel
      ),
      interpolate = TRUE,
      na.rm = TRUE
    ) +
    scale_fill_gradientn(
      colours = heat_cols,
      limits = c(0, 1),
      na.value = NA,
      guide = "none"
    )
  
  
  # Contours of focal group
  
  focal_thresholds <- thresholds[
    thresholds$Group == focal,
  ]
  
  for (level in contour_levels) {
    
    p <- p +
      geom_contour(
        data = focal_df,
        aes(
          Axis1,
          Axis2,
          z = Density,
          colour = Group
        ),
        breaks = focal_thresholds[[level]],
        linewidth = contour_width,
        show.legend = FALSE
      )
  }
  
  
  # Arrows
  
  arrows <- arrows |>
    dplyr::mutate(
      xend = Dim1 * arrow_scale,
      yend = Dim2 * arrow_scale
    )
  
  p <- p +
    geom_segment(
      data = arrows,
      aes(
        x = 0,
        y = 0,
        xend = xend,
        yend = yend
      ),
      arrow = grid::arrow(
        length = grid::unit(
          0.2,
          "cm"
        )
      ),
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = arrows,
      aes(
        x = xend * 1.05,
        y = yend * 1.05,
        label = Label
      ),
      size = 5,
      inherit.aes = FALSE
    )
  
  
  # Plot
  
  p <- p +
    scale_colour_manual(
      values = colours
    ) +
    labs(
      x = axes[1],
      y = axes[2]
    ) +
    coord_fixed(
      xlim = xlim,
      ylim = ylim,
      expand = FALSE
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      legend.position = "none",
      panel.background = element_rect(
        fill = "white",
        colour = NA
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )
  
  p
}


#### Arrows ####

numeric_arrows <- numeric_loadings |>
  dplyr::filter(
    abs(PC1) >= .5 |
      abs(PC2) >= .5
  ) |>
  dplyr::transmute(
    TraitGroup,
    Label = Trait,
    Dim1 = PC1,
    Dim2 = PC2
  )


categorical_arrows <- categorical_vectors |>
  dplyr::filter(
    abs(Cor_PCo1) >= .3 |
      abs(Cor_PCo2) >= .3
  ) |>
  dplyr::transmute(
    TraitGroup,
    Label = dplyr::recode(
      Trait,
      Phen_Start_rad = "cir_Phen_Start",
      Phen_End_rad = "cir_Phen_End"
    ),
    Dim1 = PCo1,
    Dim2 = PCo2
  )


#### Figures 1 and 2 ####

trait_groups <- c(
  AUR = 1,
  RP = 2,
  TAS = 3,
  AES = 4
)


# TPD objects

TPD_origin <- list(
  AUR = TPD_numeric$AUR$Origin,
  RP = TPD_categorical$RP$Origin,
  TAS = TPD_numeric$TAS$Origin,
  AES = TPD_categorical$AES$Origin
)

TPD_status <- list(
  AUR = TPD_numeric$AUR$Inv.Status,
  RP = TPD_categorical$RP$Inv.Status,
  TAS = TPD_numeric$TAS$Inv.Status,
  AES = TPD_categorical$AES$Inv.Status
)


# Arrows by trait group

TPD_arrows <- list(
  
  AUR = numeric_arrows |>
    dplyr::filter(
      TraitGroup == "AUR"
    ),
  
  RP = categorical_arrows |>
    dplyr::filter(
      TraitGroup == "RP"
    ),
  
  TAS = numeric_arrows |>
    dplyr::filter(
      TraitGroup == "TAS"
    ),
  
  AES = categorical_arrows |>
    dplyr::filter(
      TraitGroup == "AES"
    )
)


# Axis labels

origin_axes <- list(
  AUR = c("PC1", "PC2"),
  RP = c("PC1", "PC2"),
  TAS = c("PC1", "PC2"),
  AES = c("Dim1", "Dim2")
)

status_axes <- list(
  AUR = c("PC1", "PC2"),
  RP = c("PCo1", "PCo2"),
  TAS = c("PC1", "PC2"),
  AES = c("PCo1", "PCo2")
)


# Axis limits

x_limits <- list(
  AUR = c(-4, 4),
  RP = c(-0.35, 0.35),
  TAS = c(-4, 4),
  AES = c(-0.15, 0.15)
)

y_limits <- list(
  AUR = c(-4, 4),
  RP = c(-0.35, 0.35),
  TAS = c(-4, 4),
  AES = c(-0.15, 0.15)
)


# Arrow scaling

arrow_scale <- c(
  AUR = 3,
  RP = 0.1,
  TAS = 3,
  AES = 0.1
)

# Contour levels

origin_contour_levels <- list(
  AUR = c("t95", "t50"),
  RP = c("t95", "t50", "t20", "t05"),
  TAS = c("t95", "t50", "t20", "t05"),
  AES = c("t95", "t50", "t20", "t05")
)

status_contour_levels <- list(
  AUR = c("t95", "t50", "t20", "t05"),
  RP = c("t95", "t50", "t20", "t05"),
  TAS = c("t95", "t50", "t20", "t05"),
  AES = c("t95", "t50", "t20", "t05")
)


# Contour width

origin_contour_width <- c(
  AUR = 2,
  RP = 1.1,
  TAS = 1.1,
  AES = 1.1
)

status_contour_width <- c(
  AUR = 1.1,
  RP = 1.1,
  TAS = 1.1,
  AES = 1.1
)

for (trait_group in names(trait_groups)) {
  
  n <- trait_groups[[trait_group]]
  
  
  #Fig.1
  
  origin_plot <- patchwork::wrap_plots(
    
    tpd_panel(
      tpd = TPD_origin[[trait_group]],
      focal = "Native",
      colours = origin_cols,
      arrows = TPD_arrows[[trait_group]],
      axes = origin_axes[[trait_group]],
      xlim = x_limits[[trait_group]],
      ylim = y_limits[[trait_group]],
      arrow_scale = arrow_scale[[trait_group]],
      contour_levels = origin_contour_levels[[trait_group]],
      contour_width = origin_contour_width[[trait_group]]
    ),
    
    tpd_panel(
      tpd = TPD_origin[[trait_group]],
      focal = "Alien",
      colours = origin_cols,
      arrows = TPD_arrows[[trait_group]],
      axes = origin_axes[[trait_group]],
      xlim = x_limits[[trait_group]],
      ylim = y_limits[[trait_group]],
      arrow_scale = arrow_scale[[trait_group]],
      contour_levels = origin_contour_levels[[trait_group]],
      contour_width = origin_contour_width[[trait_group]]
    ),
    
    ncol = 2
  )
  
  
  ggsave(
    paste0(
      "Output/Fig.1/TPD_heatmap",
      n,
      "_origin.jpeg"
    ),
    origin_plot,
    width = 10,
    height = 5,
    dpi = 600,
    bg = "white"
  )
  
  
  # Fig.2
  
  status_plot <- patchwork::wrap_plots(
    
    tpd_panel(
      tpd = TPD_status[[trait_group]],
      focal = "Non-Established",
      colours = status_cols,
      arrows = TPD_arrows[[trait_group]],
      axes = status_axes[[trait_group]],
      xlim = x_limits[[trait_group]],
      ylim = y_limits[[trait_group]],
      arrow_scale = arrow_scale[[trait_group]],
      contour_levels = status_contour_levels[[trait_group]],
      contour_width = status_contour_width[[trait_group]]
    ),
    
    tpd_panel(
      tpd = TPD_status[[trait_group]],
      focal = "Naturalised",
      colours = status_cols,
      arrows = TPD_arrows[[trait_group]],
      axes = status_axes[[trait_group]],
      xlim = x_limits[[trait_group]],
      ylim = y_limits[[trait_group]],
      arrow_scale = arrow_scale[[trait_group]],
      contour_levels = status_contour_levels[[trait_group]],
      contour_width = status_contour_width[[trait_group]]
    ),
    
    tpd_panel(
      tpd = TPD_status[[trait_group]],
      focal = "Invasive",
      colours = status_cols,
      arrows = TPD_arrows[[trait_group]],
      axes = status_axes[[trait_group]],
      xlim = x_limits[[trait_group]],
      ylim = y_limits[[trait_group]],
      arrow_scale = arrow_scale[[trait_group]],
      contour_levels = status_contour_levels[[trait_group]],
      contour_width = status_contour_width[[trait_group]]
    ),
    
    ncol = 3
  )
  
  
  ggsave(
    paste0(
      "Output/Fig.2/TPD_heatmap",
      n,
      ".jpeg"
    ),
    status_plot,
    width = 10,
    height = 5,
    dpi = 600,
    bg = "white"
  )
}

#### Figure 3 ####

imp_data <- read.csv("Processed data/imp_data.csv", sep = ";", dec = ".", check.names = F)

figure3_data <- imp_data |>
  dplyr::mutate(Origin = factor(Origin, levels = c("Native", "Alien")),
                Inv.Status = factor(Inv.Status, levels = c("Native", "Non-Established",
                                                           "Naturalised", "Invasive")),
                logHeight = log(Height),
                logSLA = log(SLA)
  )

figure3_status <- figure3_data |>
  dplyr::filter(Origin == "Alien")

figure3_settings <- list(
  
  Fig3A = list(
    variable = "logSLA",
    ylab = "log(Specific leaf area)",
    letters = c(
      Native = "A",
      Alien = "B"),
    ylim = c(0,5),
    breaks = NULL
  ),
  
  Fig3B = list(
    variable = "logHeight",
    ylab = "log(Height)",
    letters = c(
      Native = "A",
      Alien = "B"
    ),
    ylim = c(-1.5, 6.5),
    breaks = NULL
  ),
  
  Fig3C = list(
    variable = "Phen_Time_End",
    ylab = "End of the flowering period (month)",
    letters = c(
      Native = "A",
      Alien = "B",
      `Non-Established` = "ab",
      Naturalised = "a",
      Invasive = "b"
    ),
    ylim = c(0, 14),
    breaks = seq(0, 12, 3)
  ),
  
  Fig3D = list(
    variable = "DS_amp",
    ylab = "Dispersal syndrome amplitude",
    letters = c(
      `Non-Established` = "a",
      Naturalised = "b",
      Invasive = "ab"
    ),
    ylim = c(.7, 7),
    breaks = 1:6
  ),
  
  Fig3E = list(
    variable = "PS_amp",
    ylab = "Pollination syndrome amplitude",
    letters = c(
      Native = "A",
      Alien = "B"
    ),
    ylim = c(.9, 2.5),
    breaks = NULL
  ),
  
  Fig3F = list(
    variable = "Tdrought",
    ylab = "Tolerance to drought",
    letters = c(
      Native = "A",
      Alien = "B",
      `Non-Established` = "a",
      Naturalised = "b",
      Invasive = "b"
    ),
    ylim = c(0, 5.7),
    breaks = 0:5
  ),
  
  Fig3G = list(
    variable = "Tshade",
    ylab = "Tolerance to shade",
    letters = c(
      `Non-Established` = "a",
      Naturalised = "b",
      Invasive = "ab"
    ),
    ylim = c(0, 5.7),
    breaks = 0:5
  ),
  
  Fig3H = list(
    variable = "LS_amp",
    ylab = "Leaf shape amplitude",
    letters = c(
      Native = "A",
      Alien = "B",
      `Non-Established` = "a",
      Naturalised = "a",
      Invasive = "b"
    ),
    ylim = c(.9, 4),
    breaks = 1:3)
  )

for (figure in names(figure3_settings)) {
  
  settings <- figure3_settings[[figure]]
  
  # Data
  
  plot_data <- dplyr::bind_rows(figure3_data |>
                                  dplyr::transmute(Group = as.character(Origin), Value = .data[[settings$variable]]),
                                figure3_data |>
                                  dplyr::filter(Origin == "Alien") |> dplyr::transmute(Group = as.character(Inv.Status), 
                                                                                       Value = .data[[settings$variable]])) |>
    dplyr::filter(!is.na(Group), is.finite(Value)) |>
    dplyr::mutate(Group = factor(Group, levels = c("Native", "Alien", "gap",
                                                   "Non-Established", "Naturalised", "Invasive")))
  
  # Boxplot
  
  plot <- ggplot(plot_data, aes(x = Group, y = Value)) + 
    geom_jitter(
    width = .15,
    height = 0,
    alpha = .20,
    size = .8,
    colour = "grey50"
    ) +
    geom_boxplot(
      outlier.shape = NA,
      fill = "white",
      colour = "grey40",
      width = .62
    ) +
    scale_x_discrete(
      drop = FALSE,
      labels = c(
        Native = "Native",
        Alien = "Alien",
        gap = "",
        `Non-Established` = "Non-Established",
        Naturalised = "Naturalised",
        Invasive = "Invasive"
      )
    ) +
    scale_y_continuous(
      breaks = if (is.null(settings$breaks)) {
        waiver()
      } else {
        settings$breaks
      },
      expand = expansion(
        mult = c(.05, .16)
      )
    ) +
    
    labs(
      x = NULL,
      y = settings$ylab
    ) +
    
    theme_bw(
      base_size = 15
    ) +
    
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(
        colour = "grey40",
        fill = NA,
        linewidth = 1.5
      ),
      axis.text.x = element_text(
        size = 12
      )
    )
  
  
  # Y-axis limits
  
  if (!is.null(settings$ylim)) {
    
    plot <- plot +
      coord_cartesian(
        ylim = settings$ylim,
        clip = "on"
      )
  }
  
  if (!is.null(settings$letters)) {
    
    y_top <- if (is.null(settings$ylim)) {
      max(
        plot_data$Value,
        na.rm = TRUE
      )
    } else {
      settings$ylim[2]
    }
    
    span <- diff(
      range(
        plot_data$Value,
        na.rm = TRUE
      )
    )
    
    letters_data <- data.frame(
      
      Group = factor(
        names(settings$letters),
        levels = levels(plot_data$Group)
      ),
      
      Label = unname(
        settings$letters
      ),
      
      y = y_top - .05 * span
    )
    
    plot <- plot +
      geom_text(
        data = letters_data,
        aes(
          x = Group,
          y = y,
          label = Label
        ),
        inherit.aes = FALSE,
        fontface = "bold",
        size = 4.5
      )
  }
  
  assign(figure, plot)
  
  ggsave(paste0("Output/Fig.3/",figure,".jpeg"),
         plot,
         width = 7.5,
         height = 5.5,
         dpi = 600,
         bg = "white")
}
