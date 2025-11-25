#Updated Script as of 11/20/25
#need micro_storm_data.csv to be located within the same project, in a folder called 'data'

####Set Up####
rm(list=ls())

#Load packages
library("readr")
library("devtools")
library("factoextra")
library("dplyr")
library("tidyr")
library("vegan")
library("googlesheets4")
library(car)
library(ggplot2)
library("forcats")
library("lmPerm")
library(ggrepel)
library(tibble)
library(patchwork)

#reading data from .csv file
micro_data <- read_csv("data/micro_storm_data.csv") #micro_data = full storm dataset

#Rename Labels
micro_data<- micro_data %>% 
  mutate (treatbysite = paste(treatment, microhabitat)) %>% #creating new column with site and pre/post storm ID
  mutate(treatbysite = factor(
    case_when(
      treatbysite == "post_storm reefcrest" ~ "Post-Storm Reef Crest",
      treatbysite == "post_storm reefflat" ~ "Post-Storm Reef Flat",
      treatbysite == "post_storm reefslope" ~ "Post-Storm Reef Base", #changing slope to base
      treatbysite == "post_storm sandflat" ~ "Post-Storm Sand Flat",
      treatbysite == "pre_storm reefcrest" ~ "Pre-Storm Reef Crest",
      treatbysite == "pre_storm reefflat" ~ "Pre-Storm Reef Flat",
      treatbysite == "pre_storm reefslope" ~ "Pre-Storm Reef Base", #changing slope to base
      treatbysite == "pre_storm sandflat" ~ "Pre-Storm Sand Flat",
      TRUE ~ treatbysite
    ),
    levels = c( #ordering factors
      "Pre-Storm Reef Crest",
      "Pre-Storm Reef Flat",
      "Pre-Storm Reef Base", 
      "Pre-Storm Sand Flat",
      "Post-Storm Reef Crest",
      "Post-Storm Reef Flat",
      "Post-Storm Reef Base",
      "Post-Storm Sand Flat"
    )
  )
  ) %>%
  rename( #editing trait labels
    "SA:DW" = sa_dw,
    "H" = h,
    "H:V" = h_v,
    "H:WW" = h_ww,
    "SA:V" = sa_v, 
    "TS" = ts,
    "T" = t,
    "% Calc" = p_calcification,
    "DW:WW" = dw_ww
  )

####Full Storm PERMANOVA####
# Create trait matrix
traits.mat <- as.matrix(micro_data[, 5:13]) 

# Calculate Bray-Curtis dissimilarity
traits.dist <- vegdist(traits.mat, method = "bray")

# Run PERMANOVA with main effects and interaction
set.seed(36)
traits.div <- adonis2(traits.dist ~ microhabitat * treatment,
                      data = micro_data,
                      permutations = 999,
                      method = "bray")

# View results
traits.div

####Full Storm PCA####
#making PCA for full storm dataset
microstormPCA <- prcomp(micro_data[, 5:13],  scale = TRUE) #5-13 is all of the traits so only want those columns
summary(microstormPCA)

pca_scores <- as.data.frame(microstormPCA$x)
pca_scores$Group <- micro_data$treatbysite  # Grouped by Treatment:Habitat Key

# Create PCA scores dataframe with IDs and Groups (using Key column)
pca_scores_with_id <- bind_cols(
  micro_data %>% select(id, treatbysite),  # Keep the ID and Key columns for reference
  as.data.frame(microstormPCA$x)   # Add the PCA scores
)

# Eigenvalues
eigenvalues <- microstormPCA$sdev^2
# Assign names to each eigenvalue
names(eigenvalues) <- paste0("PC", seq_along(eigenvalues))
eigenvalues

# # Plot basic PCA with shaded ellipses
# ggplot(pca_scores_with_id, aes(x = PC1, y = PC2, color = Key, fill = Key)) +
#   geom_point(size = 3) +  # Points for each sample
#   stat_ellipse(aes(color = Key), level = 0.95,
#                geom = "polygon",  # Shading the ellipses
#                alpha = 0.3) +  # Transparency for ellipses
#   labs(title = "Basic PCA Visualization", x = "PC1", y = "PC2") +
#   theme_minimal()

# Extract loadings (PC1 and PC2)
loadings <- as.data.frame(microstormPCA$rotation[, 1:2])
loadings$trait <- rownames(loadings)

# Compute label angles
loadings$angle <- atan2(loadings$PC2, loadings$PC1) * 180 / pi
loadings$angle <- ifelse(loadings$angle > 90 | loadings$angle < -90,
                         loadings$angle + 180,
                         loadings$angle)



# Arrow scale (same as before)
arrow_scale <- 6

# Arrow endpoints
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale

#scaling labels to make them a good distance from arrow
loadings$label_scale <- 1  
loadings$label_scale[loadings$trait == "H:WW"] <- 1.25
loadings$label_scale[loadings$trait == "H:V"]  <- 1.15
loadings$label_scale[loadings$trait == "% Calc"] <- 1.35
loadings$label_scale[loadings$trait == "DW:WW"] <- 1.35
loadings$label_scale[loadings$trait == "SA:V"] <- 1.35
loadings$label_scale[loadings$trait == "T"] <- 1.15
loadings$label_scale[loadings$trait == "SA:DW"] <- 1.25
loadings$label_scale[loadings$trait == "TS"] <- 1.25
loadings$label_scale[loadings$trait == "H"] <- 1.1

# Label position = arrow endpoint × label_scale
loadings$label_x <- loadings$xend * loadings$label_scale
loadings$label_y <- loadings$yend * loadings$label_scale


# Manually adjusting labels for 'SA:V' and 'T' so they are not stacked on top of each other (spacing them out)
loadings$label_shift <- 0  # Default label shift (no shift)
# Add manual shifts for 'SA:V' and 'T' because they are too close together
loadings$label_y[loadings$trait == "SA:V"] <- loadings$label_y[loadings$trait == "SA:V"] + 0.2
loadings$label_y[loadings$trait == "T"]    <- loadings$label_y[loadings$trait == "T"] - 0.2

#Visualizing PCA Plot (template - this figure is not actually used in paper, just for reference)
ggplot(pca_scores_with_id, aes(x = PC1, y = PC2, color = treatbysite, fill = treatbysite, label=id)) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.5)+
  # geom_point(size = 1.25) +  # Points for each sample
  geom_text_repel(size = 3, show.legend = FALSE) +
  stat_ellipse(aes(color = treatbysite), level = 0.5, #50% confidence interval 
               geom = "polygon",  # Shading the ellipses
               alpha = 0.3) +  # Transparency for ellipses
  # Arrows from origin
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    arrow = arrow(type = "closed", length = unit(0.075, "inches")),
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Trait labels perfectly aligned with arrow direction
  geom_text(
    data = loadings,
    aes(x = label_x, y = label_y, label = trait, angle = angle),
    color = "black",
    size = 4,
    hjust = 0.5,
    vjust = 0.5,
    inherit.aes = FALSE
  ) +
  labs(title = "PCA with Trait Contribution Arrows", x = "PC1 (26.2%)", y = "PC2 (19.2%)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  scale_color_manual(values = c("Pre-Storm Sand Flat" = "#d8bdea",
                                "Pre-Storm Reef Flat" = "#c1ea76",
                                "Pre-Storm Reef Crest" = "#f3c8c3",
                                "Pre-Storm Reef Base" = "#a1e5e8",
                                "Post-Storm Sand Flat" = "#ba57ff",
                                "Post-Storm Reef Flat" = "#69a400",
                                "Post-Storm Reef Crest" = "#f55d56",
                                "Post-Storm Reef Base" = "#00b3b8")) +  # Custom color palette
  scale_fill_manual(values = c("Pre-Storm Sand Flat" = "#d8bdea",
                               "Pre-Storm Reef Flat" = "#c1ea76",
                               "Pre-Storm Reef Crest" = "#f3c8c3",
                               "Pre-Storm Reef Base" = "#a1e5e8",
                               "Post-Storm Sand Flat" = "#ba57ff",
                               "Post-Storm Reef Flat" = "#69a400",
                               "Post-Storm Reef Crest" = "#f55d56",
                               "Post-Storm Reef Base" = "#00b3b8")) 


####PCA Plot Function####
create_pca_plot <- function(data, highlight_groups, color_values, title) {
  ggplot(pca_scores_with_id, aes(x = PC1, y = PC2, color = treatbysite, fill = treatbysite)) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
    geom_vline(xintercept = 0, color = "gray70", linewidth = 0.5)+
    geom_point(size = 1.25) +  # Points for each sample
    stat_ellipse(data = data %>% filter(treatbysite %in% highlight_groups), 
                 aes(color = treatbysite), level = 0.5, #50% confidence interval 
                 geom = "polygon",  # Shading the ellipses
                 alpha = 0.3) +  # Transparency for ellipses
    # Arrows from origin
    geom_segment(
      data = loadings,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      arrow = arrow(type = "closed", length = unit(0.075, "inches")),
      color = "black",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    
    # Trait labels perfectly aligned with arrow direction
    geom_text(
      data = loadings,
      aes(x = label_x, y = label_y, label = trait, angle = angle),
      color = "black",
      size = 4,
      hjust = 0.5,
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    coord_fixed()+
        labs(title = title, x = "PC1 (26.2%)", y = "PC2 (19.2%)") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank()   # Remove minor gridlines
    ) +
    scale_color_manual(values = c("Pre-Storm Sand Flat" = "#d8bdea",
                                  "Pre-Storm Reef Flat" = "#c1ea76",
                                  "Pre-Storm Reef Crest" = "#f3c8c3",
                                  "Pre-Storm Reef Base" = "#a1e5e8",
                                  "Post-Storm Sand Flat" = "#ba57ff",
                                  "Post-Storm Reef Flat" = "#69a400",
                                  "Post-Storm Reef Crest" = "#f55d56",
                                  "Post-Storm Reef Base" = "#00b3b8")) +  # Custom color palette
    scale_fill_manual(values = c("Pre-Storm Sand Flat" = "#d8bdea",
                                 "Pre-Storm Reef Flat" = "#c1ea76",
                                 "Pre-Storm Reef Crest" = "#f3c8c3",
                                 "Pre-Storm Reef Base" = "#a1e5e8",
                                 "Post-Storm Sand Flat" = "#ba57ff",
                                 "Post-Storm Reef Flat" = "#69a400",
                                 "Post-Storm Reef Crest" = "#f55d56",
                                 "Post-Storm Reef Base" = "#00b3b8"))  +
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    theme(legend.position = "none") #removes legend if needed
}

####Microhabitat PCA Plots####
#Reef Crest
highlight_rc <- c("Pre-Storm Reef Crest", "Post-Storm Reef Crest")
rcpca<- create_pca_plot(pca_scores_with_id, highlight_rc, color_values, 
                        NULL)
rcpca
#Reef Flat
highlight_rf <- c("Pre-Storm Reef Flat", "Post-Storm Reef Flat")
rfpca<- create_pca_plot(pca_scores_with_id, highlight_rf, color_values, 
                        NULL)
rfpca

#Reef Base
highlight_rb <- c("Pre-Storm Reef Base", "Post-Storm Reef Base")
rbpca<-create_pca_plot(pca_scores_with_id, highlight_rb, color_values, 
                       NULL)

#Sand Flat
highlight_sf <- c("Pre-Storm Sand Flat", "Post-Storm Sand Flat")
sfpca<- create_pca_plot(pca_scores_with_id, highlight_sf, color_values, 
                        NULL)

pca_gridded <- (rcpca | rfpca) /
  (rbpca | sfpca)
pca_gridded

####% Trait Contributions for Full Storm####
loadings_abs <- loadings %>%
  select(trait, PC1, PC2) %>%
  mutate(across(c(PC1, PC2), abs)) 

# Square the loadings
loadings_squared <- loadings %>%
  select(trait, PC1, PC2) %>%
  mutate(across(c(PC1, PC2), ~ .^2))

# Calculate percent contribution
loadings_percent <- loadings_squared %>%
  mutate(PC1 = 100 * PC1 / sum(PC1),
         PC2 = 100 * PC2 / sum(PC2))

# Sort for plotting
pc1_sorted <- loadings_percent[order(-loadings_percent$PC1), c("trait", "PC1")]
pc2_sorted <- loadings_percent[order(-loadings_percent$PC2), c("trait", "PC2")]

# Calculate null model line
n_traits <- nrow(loadings_percent)
null_contribution <- 100 / n_traits

# % Contribution Barplots with null line
par(mfrow = c(1, 2), mar = c(7, 4, 4, 2))  # Increased bottom margin
global_max <- max(c(pc1_sorted$PC1, pc2_sorted$PC2))
ylim_shared <- c(0, global_max * 1.1)   # add 10% headroom
barplot(pc1_sorted$PC1,
        names.arg = pc1_sorted$trait,
        col = "skyblue",
        las = 2,
        xlab = "",
        ylim = ylim_shared)
mtext("Trait", side = 1, line = 4) 
mtext("PC1 (% contribution)", side = 2, line = 2.3) 
abline(h = null_contribution, lty = 2, col = "red")

barplot(pc2_sorted$PC2,
        names.arg = pc2_sorted$trait,
        col = "skyblue",
        las = 2,
        xlab = "",
        ylim = ylim_shared)
mtext("Trait", side = 1, line = 4) 
mtext("PC2 (% contribution)", side = 2, line = 2.3) 
abline(h = null_contribution, lty = 2, col = "red")

####Pre Storm only PERMANOVA####
#One Factor Permanova
prestorm_data <- micro_data %>%
  filter(treatment == "pre_storm")

# Create matrix directly
traits.mat <- as.matrix(prestorm_data[, 5:13]) 

# Calculate Bray-Curtis dissimilarity
traits.dist <- vegdist(traits.mat, method = 'bray')

# Run PERMANOVA
set.seed(36)
traits.div <- adonis2(traits.dist ~ microhabitat,
                      data = prestorm_data,
                      permutations = 999,
                      method = "bray")
# print results
traits.div

####Just Pre-Storm Data PCA####
prestormPCA <- prcomp(prestorm_data[, 5:13], scale = TRUE)
prestorm_scores <- as.data.frame(prestormPCA$x)
prestorm_scores$Group <- prestorm_data$treatbysite  # Grouped by Treatment:Habitat Key

# PCA scores (with ID + treatbysite)
pre_scores_id <- bind_cols(
  prestorm_data %>% select(id, treatbysite),
  as.data.frame(prestormPCA$x)
)

# Loadings
loadings_pre <- as.data.frame(prestormPCA$rotation[, 1:2])
loadings_pre$trait <- rownames(loadings_pre)

# Arrow angle
loadings_pre$angle <- atan2(loadings_pre$PC2, loadings_pre$PC1) * 180 / pi

# Flip upside-down labels
loadings_pre$angle <- ifelse(
  loadings_pre$angle > 90 | loadings_pre$angle < -90,
  loadings_pre$angle + 180,
  loadings_pre$angle
)

#adding on arrows
arrow_scale_pre <- 5
# Arrow endpoints
loadings_pre$xend <- loadings_pre$PC1 * arrow_scale_pre
loadings_pre$yend <- loadings_pre$PC2 * arrow_scale_pre

#scaling labels to make them a good distance from arrow
loadings_pre$label_scale <- 1  
loadings_pre$label_scale[loadings_pre$trait == "H:WW"] <- 1.15
loadings_pre$label_scale[loadings_pre$trait == "H:V"]  <- 1.15
loadings_pre$label_scale[loadings_pre$trait == "% Calc"] <- 1.35
loadings_pre$label_scale[loadings_pre$trait == "DW:WW"] <- 1.3
loadings_pre$label_scale[loadings_pre$trait == "SA:V"] <- 1.25
loadings_pre$label_scale[loadings_pre$trait == "T"] <- 1.15
loadings_pre$label_scale[loadings_pre$trait == "SA:DW"] <- 1.2
loadings_pre$label_scale[loadings_pre$trait == "TS"] <- 1.75
loadings_pre$label_scale[loadings_pre$trait == "H"] <- 1.1

# Label position = arrow endpoint × label_scale
loadings_pre$label_x <- loadings_pre$xend * loadings_pre$label_scale
loadings_pre$label_y <- loadings_pre$yend * loadings_pre$label_scale

# Calculate mean trait value for each sample
centroids_pre <- pre_scores_id %>%
  group_by(treatbysite) %>%
  summarize(
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2)
  )

#prestorm pca plot
ggplot(pre_scores_id, aes(x = PC1, y = PC2, color = treatbysite, fill = treatbysite)) +
  geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray70", linewidth = 0.5)+
  geom_point() +
  geom_point(
    data = centroids_pre,
    aes(x = mean_PC1, y = mean_PC2, color = treatbysite),
    size = 4,        # large
    shape = 16,      # solid circle
    inherit.aes = FALSE
  )+
  stat_ellipse(aes(color = treatbysite), level = 0.50, #50% confidence interval 
               geom = "polygon",  # Shading the ellipses
               alpha = 0.3) +  # Transparency for ellipses
  # Arrows from origin showing trait contribution (override aesthetics)
  geom_segment(data = loadings_pre,
               aes(x = 0, y = 0, xend = PC1 * arrow_scale_pre, yend = PC2 * arrow_scale_pre),
               arrow = arrow(type = "closed", length = unit(0.075, "inches")),
               color = "black", linewidth = .5, inherit.aes = FALSE) +
  # Trait labels at end of arrows
  geom_text(
    data = loadings_pre,
    aes(x = label_x, y = label_y, label = trait, angle = angle),
    color = "black",
    hjust = 0.5,
    vjust = 0.5,
    size = 4,
    inherit.aes = FALSE
  ) +
  labs(title = NULL, x = "PC1 (23.9%)", y = "PC2 (21.1%)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank()   # Remove minor gridlines
  ) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.1)) +
  scale_color_manual(name = "Key", values = c("Pre-Storm Sand Flat" = "#ba57ff", #made all prestorm outlines dark for the prestorm only
                                "Pre-Storm Reef Flat" = "#69a400",
                                "Pre-Storm Reef Crest" = "#f55d56",
                                "Pre-Storm Reef Base" = "#00b3b8",
                                "Post-Storm Sand Flat" = "#ba57ff",
                                "Post-Storm Reef Flat" = "#69a400",
                                "Post-Storm Reef Crest" = "#f55d56",
                                "Post-Storm Reef Base" = "#00b3b8")) +  # Custom color palette
  scale_fill_manual(name = "Key", values = c("Pre-Storm Sand Flat" = "#d8bdea",
                               "Pre-Storm Reef Flat" = "#c1ea76",
                               "Pre-Storm Reef Crest" = "#f3c8c3",
                               "Pre-Storm Reef Base" = "#a1e5e8",
                               "Post-Storm Sand Flat" = "#ba57ff",
                               "Post-Storm Reef Flat" = "#69a400",
                               "Post-Storm Reef Crest" = "#f55d56",
                               "Post-Storm Reef Base" = "#00b3b8"))  

####% Contributions for Just Pre- Storm####
loadings_pre_abs <- loadings_pre %>%
  select(trait, PC1, PC2) %>%
  mutate(across(c(PC1, PC2), abs)) 

# Square the loadings
loadings_pre_squared <- loadings_pre %>%
  select(trait, PC1, PC2) %>%
  mutate(across(c(PC1, PC2), ~ .^2))

# Calculate percent contribution
loadings_pre_percent <- loadings_pre_squared %>%
  mutate(PC1 = 100 * PC1 / sum(PC1),
         PC2 = 100 * PC2 / sum(PC2))

# Sort for plotting
pc1_pre_sorted <- loadings_pre_percent[order(-loadings_pre_percent$PC1), c("trait", "PC1")]
pc2_pre_sorted <- loadings_pre_percent[order(-loadings_pre_percent$PC2), c("trait", "PC2")]

# Calculate null model line
n_traits <- nrow(loadings_pre_percent)
null_contribution <- 100 / n_traits

# % Contribution Barplots with null line
global_max_pre <- max(c(pc1_pre_sorted$PC1, pc2_pre_sorted$PC2))
ylim_shared <- c(0, global_max_pre * 1.1)   # add 10% headroom
par(mfrow = c(1, 2), mar = c(7, 4, 4, 2))  # Increased bottom margin

# For two plots side by side
barplot(pc1_pre_sorted$PC1,
        names.arg = pc1_pre_sorted$trait,
        col = "skyblue",
        las = 2,
        ylim = ylim_shared,
        xlab = "")  # Suppress default x-axis label
mtext("Trait", side = 1, line = 4.5) 
mtext("PC1 (% contribution)", side = 2, line = 2.3) 
abline(h = null_contribution, lty = 2, col = "red")

barplot(pc2_pre_sorted$PC2,
        names.arg = pc2_pre_sorted$trait,
        col = "skyblue",
        las = 2,
        ylim = ylim_shared,
        xlab = "")
mtext("Trait", side = 1, line = 4.5) 
mtext("PC2 (% contribution)", side = 2, line = 2.3) 
abline(h = null_contribution, lty = 2, col = "red")

####Bar Graphs####
#clean up titles/reorder Factors
micro_data$microhabitat_f <- factor(micro_data$microhabitat,
                                    levels = c("sandflat", "reefflat", "reefcrest", "reefslope"),
                                    labels = c("Sand Flat", "Reef Flat", "Reef Crest", "Reef Base"))
#Function to make bar graphs
micro_barplot <- function(data, yvar, color1, color2, yylab, show_legend = TRUE) {
  
  # Summarize: mean, SD, count, standard error
  micro_summary_stats <- data %>%
    group_by(microhabitat, treatment) %>%
    summarize(across(all_of(yvar), list(mean = mean, sd = sd, count = length), .names = "{col}_{fn}"),
              .groups = "drop") %>%
    mutate(
      mean = !!sym(paste0(yvar, "_mean")),
      sd = !!sym(paste0(yvar, "_sd")),
      count = !!sym(paste0(yvar, "_count")),
      st_error = sd / sqrt(count)
    )
  
  # Plot
  Micro_Plot <- ggplot(micro_summary_stats,
                       aes(x = factor(microhabitat, levels = c('sandflat', 'reefflat', 'reefcrest', 'reefslope')),
                           y = mean,
                           fill = fct_rev(treatment),
                           color = fct_rev(treatment))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) + 
    geom_errorbar(aes(ymin = mean - st_error, ymax = mean + st_error),
                  color = "black", position = position_dodge(width = 0.8), width = 0.2) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = ifelse(show_legend, "right", "none")) +
    xlab("Microhabitat") + 
    ylab(yylab) +
    scale_fill_manual(values = c(color1, color2), labels = c("Pre Storm", "Post Storm")) + 
    scale_color_manual(values = c(color1, color2), labels = c("Pre Storm", "Post Storm")) + 
    scale_x_discrete(labels = c("Sand Flat", "Reef Flat", "Reef Crest", "Reef Slope"))
  
  return(Micro_Plot)
}

#Bar Graphs for Each Trait
tough_bar <- micro_barplot(data=micro_data, yvar = 'T', color2="#1B9E77", color1 = "#56B4E9", 
                         yylab = "Toughness", show_legend = FALSE)
tough_bar

p_calc_bar <- micro_barplot(data=micro_data, yvar = '% Calc', color2="#1B9E77", color1 = "#56B4E9", 
                       yylab = "Percent Calcification (%)", show_legend = FALSE)
p_calc_bar

dw_ww_bar <- micro_barplot(data=micro_data, yvar = 'DW:WW', color2="#1B9E77", color1 = "#56B4E9", 
                           yylab = "Dry Weight to Wet Weight", show_legend = FALSE)
dw_ww_bar

tensile_bar <- micro_barplot(data=micro_data, yvar = 'TS', color2="#1B9E77", color1 = "#56B4E9", 
                          yylab = "Tensile Strength", show_legend = FALSE)
tensile_bar

height_bar <- micro_barplot(data=micro_data, yvar = 'H', color2="#1B9E77", color1 = "#56B4E9", 
                             yylab = "Height", show_legend = FALSE)
height_bar

hv_bar <- micro_barplot(data=micro_data, yvar = 'H:V', color2="#1B9E77", color1 = "#56B4E9", 
                            yylab = "Height to Volume", show_legend = FALSE)
hv_bar

h_ww_bar <- micro_barplot(data=micro_data, yvar = 'H:WW', color2="#1B9E77", color1 = "#56B4E9", 
                            yylab = "Height to Wet Weight", show_legend = FALSE)
h_ww_bar

sa_v_bar <- micro_barplot(data=micro_data, yvar = 'SA:V', color2="#1B9E77", color1 = "#56B4E9", 
                          yylab = "Surface Area to Volume", show_legend = FALSE)
sa_v_bar

sa_dw_bar <- micro_barplot(data=micro_data, yvar = 'SA:DW', color2="#1B9E77", color1 = "#56B4E9", 
                          yylab = "Surface Area to Dry Weight", show_legend = FALSE)
sa_dw_bar

####Univariate Trait Analysis####

#Toughness
tough_aov <- aov(T ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(tough_aov)) #tests normality
leveneTest(T ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#% Calcification
p_calc_aov <- aov(`% Calc` ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(p_calc_aov)) #tests normality
leveneTest(`% Calc` ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance
summary(p_calc_aov)

#DW:WW
dw_ww_aov <- aov(`DW:WW` ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(dw_ww_aov)) #tests normality
leveneTest(`DW:WW` ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#Tensile Strength
ts_aov <- aov(TS ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(ts_aov)) #tests normality
leveneTest(TS ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#Height
h_aov <- aov(H ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(h_aov)) #tests normality
leveneTest(H ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#H:V
h_v_aov <- aov(`H:V` ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(h_v_aov)) #tests normality
leveneTest(`H:V` ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#H:WW
h_ww_aov <- aov(`H:WW` ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(h_ww_aov)) #tests normality
leveneTest(`H:WW` ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#SA:V
sa_v_aov <- aov(`SA:V` ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(sa_v_aov)) #tests normality
leveneTest(`SA:V` ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance

#SA:DW
sa_dw_aov <- aov(`SA:DW` ~ microhabitat * treatment, data = micro_data)
shapiro.test(residuals(sa_dw_aov)) #tests normality
leveneTest(`SA:DW` ~ microhabitat * treatment, data = micro_data) #tests homogeneity of variance


