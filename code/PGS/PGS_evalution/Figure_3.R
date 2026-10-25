# Setup #
library(pacman)
p_load(tidyverse)


###### UKB ############

root = "/lindstroem/austin_working/Dissertation/Aim_3/Rare_variant_PGS/UKB_trained/"
e_folder = "UKB_external/"

# Load data #
e_valid_breast = read.csv(paste0(root, e_folder, "UKB_trained_breast_external_OR_by_quantile_table.txt")) %>%
  mutate(Cancer = "Breast cancer", Validation = "AoU (External)")
e_valid_colorectal = read.csv(paste0(root, e_folder, "UKB_trained_colorectal_external_OR_by_quantile_table.txt"))  %>%
  mutate(Cancer = "Colorectal cancer", Validation = "AoU (External)")
e_valid_prostate = read.csv(paste0(root, e_folder, "UKB_trained_prostate_external_OR_by_quantile_table.txt"))  %>%
  mutate(Cancer = "Prostate cancer", Validation = "AoU (External)")

validation = rbind.data.frame( e_valid_breast, e_valid_colorectal, e_valid_prostate) %>%
  filter(Superpopulation == "EUR")

validation$Cancer = factor(validation$Cancer, levels = c("Breast cancer", "Colorectal cancer", "Prostate cancer"))

# Plot #
validation <- validation %>%
  mutate(
    # CV quantile (before comma)
    cv_lab = sub(",.*$", "", JointLabel),
    
    # RV quantile (after comma)
    rv_lab = trimws(sub("^[^,]*,\\s*", "", JointLabel)),
    rv_lab = gsub("\\[90-100%\\)", "[90-100%]", rv_lab),
    
    # Reference ONLY when:
    #  - RV lower bound is 0
    #  - CV quantile is [33-66%)
    is_ref = grepl("^\\[0", rv_lab) & cv_lab == "[33–66%)",
    
    # X-axis label
    x_lab = cv_lab,
    
    # Color label with special handling
    color_lab = dplyr::case_when(
      is_ref ~ "Reference (non-cRV carriers)",
      grepl("^\\[0", rv_lab) ~ "non-cRV carriers",
      grepl("75%\\)$", rv_lab) ~ "cRV carriers < 75%",
      TRUE ~ rv_lab
    )
  )


rv_levels <- validation %>%
  filter(!is_ref) %>%
  distinct(color_lab, rv_lab) %>%
  mutate(
    lower = dplyr::case_when(
      color_lab == "non-carriers" ~ -Inf,
      TRUE ~ as.numeric(sub("^\\[(\\d+).*", "\\1", rv_lab))
    )
  ) %>%
  arrange(lower) %>%
  pull(color_lab)



rv_color_map <- setNames(
  c(
    "#084594",
    "#084594",
    "#084594",
    "#00441B",
    "#00441B",
    "#00441B",
    "darkorange",
    "red4"
  )[seq_along(rv_levels)],
  rv_levels
)


# Add reference color
rv_color_map["Reference (non-cRV carriers)"] <- "black"

# Ensure unique ordered levels (Reference first conceptually)
ordered_levels <- unique(c("Reference (non-cRV carriers)", rv_levels))

validation <- validation %>%
  mutate(
    color_lab = factor(color_lab, levels = ordered_levels)
  ) %>%
  group_by(x_lab) %>%
  mutate(
    plot_index = dplyr::case_when(
      is_ref ~ 1L,
      TRUE ~ NA_integer_
    )
  ) %>%
  mutate(
    plot_index = ifelse(
      is.na(plot_index),
      as.integer(color_lab) + 1L,
      plot_index
    )
  ) %>%
  ungroup()


# Plot
plot <- ggplot(
  validation,
  aes(
    x = x_lab,
    y = OR,
    ymin = OR_lower,
    ymax = OR_upper,
    color = color_lab,
    group = interaction(Cancer, x_lab, plot_index, Validation)
  )
) +
  geom_pointrange(
    position = position_dodge(width = 0.7),
    size = 1.2
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray30",
    linewidth = 0.8
  ) +
  labs(
    x = "CV PGS quantile",
    y = "OR",
    color = "RV PGS quantile"
  ) +
  scale_y_log10(
    breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 10),
    minor_breaks = NULL
  ) +
  scale_color_manual(values = rv_color_map, drop = FALSE) +
  scale_shape_manual(values = c(15, 16)) +
  theme_classic(base_size = 20) +
  theme(
    # Axis titles
    axis.title = element_text(size = 22),
    
    # Axis tick labels
    axis.text.x = element_text(
      size = 17,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 18),
    
    # Legend
    legend.position = "top",
    legend.title = element_text(size = 19),
    legend.text = element_text(size = 19),
    legend.key.size = unit(1.2, "cm"),
    
    # Facet labels
    strip.background = element_blank(),
    strip.text = element_text(size = 20, face = "bold")
  ) +
  facet_wrap(
    ~ Cancer,
    nrow = 3,
    ncol = 1,
    scales = "free_y"
  )


ggsave("/lindstroem/austin_working/RV_h2_and_PGS_manuscript_revision/Manuscript_Figure_3.png", plot, width = 18, height = 16, dpi = 800)
