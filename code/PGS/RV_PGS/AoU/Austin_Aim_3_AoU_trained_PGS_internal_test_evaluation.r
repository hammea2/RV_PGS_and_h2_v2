#install.packages("pacman")
library(pacman)

p_load(ggplot2)

p_load(patchwork)

p_load(rlang, dplyr, data.table,
       reactable, R.utils,
       purrr, pROC, DescTools,
       prettydoc, furrr, broom,
       forcats,
       knitr, kableExtra)


breast_data = fread("./Aim_3/Rare_variant_PGS/RICE/AoU_RV_breast_PGS_internal_testing_dataset.txt")
colorectal_data = fread("./Aim_3/Rare_variant_PGS/RICE/AoU_RV_colorectal_PGS_internal_testing_dataset.txt")
prostate_data = fread("./Aim_3/Rare_variant_PGS/RICE/AoU_RV_prostate_PGS_internal_testing_dataset.txt")

breast_CV_PGS_score = fread("./Aim_3/Common_variant_PGS/scores/breast_scores.sscore")
colorectal_CV_PGS_score = fread("./Aim_3/Common_variant_PGS/scores/colorectal_scores.sscore")
prostate_CV_PGS_score = fread("./Aim_3/Common_variant_PGS/scores/prostate_scores.sscore")

breast_raw_data = merge(breast_data, breast_CV_PGS_score, by.x = "person_id", by.y = "IID") %>%
  select(names(breast_data), SCORE1_AVG)
colorectal_raw_data = merge(colorectal_data, colorectal_CV_PGS_score, by.x = "person_id", by.y = "IID") %>%
  select(names(colorectal_data), SCORE1_AVG)
prostate_raw_data = merge(prostate_data, prostate_CV_PGS_score, by.x = "person_id", by.y = "IID") %>%
  select(names(prostate_data), SCORE1_AVG)

standardize_pgs_pca <- function(data, pgs_col, pca_cols, ref_filter = NULL) {
  # --- 1. Define reference population ---
  if (!is.null(ref_filter)) {
    ref_data <- data %>% filter(!!parse_expr(ref_filter))
  } else {
    ref_data <- data
  }
  
  # --- 2. Fit regression model for mean adjustment (Z_norm1) ---
  formula_mean <- as.formula(paste(pgs_col, "~", paste(pca_cols, collapse = " + ")))
  mean_model <- lm(formula_mean, data = ref_data)
  
  # Predicted mean PGS for all individuals
  data$pred_mean <- predict(mean_model, newdata = data)
  
  # Residuals
  data$resid_mean <- data[[pgs_col]] - data$pred_mean
  
  # SD of residuals in reference population
  if (!is.null(ref_filter)) {
    sd_ref <- sd(data$resid_mean[with(data, eval(parse(text = ref_filter)))], na.rm = TRUE)
  } else {
   sd_ref <- sd(data$resid_mean, na.rm = TRUE)
  }

  # Z_norm1: mean-adjusted
  data$Z_norm1 <- data$resid_mean / sd_ref
  
  # --- 3. Fit regression model for variance adjustment (Z_norm2) ---
  # Create ref_data with residuals for variance model
  ref_data <- if (!is.null(ref_filter)) {
    data %>% filter(!!parse_expr(ref_filter))
  } else {
  data
  }
  
  ref_data <- ref_data %>% mutate(resid_sq = resid_mean^2)
  formula_sd <- as.formula(paste("resid_sq ~", paste(pca_cols, collapse = " + ")))
  sd_model <- lm(formula_sd, data = ref_data)
  
  # Predicted SD for all individuals
  data$pred_var <- predict(sd_model, newdata = data)
  data$pred_sd <- sqrt(pmax(data$pred_var, 1e-8))
  
  # Z_norm2: mean + variance normalized
  data$Z_norm2 <- data$resid_mean / data$pred_sd
  
  return(data)
}

evaluate_auc_r2_two_scores <- function(
  data,
  score1_col,
  score2_col,
  covar_col = NULL, 
  case_col,
  ancestry_col,
  pops = NULL,
  ci_level = 0.95
) {
  
  # Filter populations if requested
  if (!is.null(pops)) {
    data <- data %>% filter(.data[[ancestry_col]] %in% pops)
  }
  
  # Function to fit a model and extract AUC + R²
  fit_and_extract <- function(sub, formula_str, model_label, resp, ci_level, score_names) {
    
    model <- glm(
      as.formula(formula_str),
      data = sub,
      family = binomial(link = "logit")
    )
    
    lp <- predict(model, type = "link")
    
    roc_obj <- tryCatch(
      roc(response = resp, predictor = lp),
      error = function(e) NULL
    )
    
    if (!is.null(roc_obj)) {
      ci <- ci.auc(roc_obj, conf.level = ci_level)
      auc       <- round(as.numeric(roc_obj$auc), 3)
      auc_lower <- round(ci[1], 3)
      auc_upper <- round(ci[3], 3)
    } else {
      auc <- auc_lower <- auc_upper <- NA_real_
    }
    
    pseudo_r2 <- round(PseudoR2(model, which = "Nagelkerke"), 5)
    
    tibble(
      Model = model_label,
      Score = score_names,
      auc = auc,
      auc_lower = auc_lower,
      auc_upper = auc_upper,
      nagelkerke_R2 = pseudo_r2
    )
  }
  
  # Process by ancestry group
  results <- data %>%
    group_split(across(all_of(ancestry_col))) %>%
    map_dfr(function(sub) {
      ancestry <- unique(sub[[ancestry_col]])
      if (nrow(sub) == 0) return(NULL)
      
      resp <- sub[[case_col]]
      
      # skip if only one outcome class
      if (length(unique(resp)) <= 1) {
        return(tibble(
          Superpopulation = ancestry,
          n = nrow(sub),
          n_cases = sum(resp == 1, na.rm = TRUE),
          n_controls = sum(resp == 0, na.rm = TRUE),
          Model = c("Combined","Score1_only","Score2_only"),
          Score = c(paste(score1_col, "+", score2_col), score1_col, score2_col),
          auc = NA_real_,
          auc_lower = NA_real_,
          auc_upper = NA_real_,
          nagelkerke_R2 = NA_real_
        ))
      }
      
      #‒‒‒‒ FIT ALL THREE MODELS ‒‒‒‒
      combined <- fit_and_extract(
        sub,
        paste(case_col, "~", score1_col, "+", score2_col),
        "Combined",
        resp,
        ci_level,
        score_names = paste(score1_col, "+", score2_col)
      )
      
      score1 <- fit_and_extract(
        sub,
        paste(case_col, "~", score1_col),
        "Score1_only",
        resp,
        ci_level,
        score_names = score1_col
      )
      
      score2 <- fit_and_extract(
        sub,
        paste(case_col, "~", score2_col),
        "Score2_only",
        resp,
        ci_level,
        score_names = score2_col
      )
      
      bind_rows(combined, score1, score2) %>%
        mutate(
          Superpopulation = ancestry,
          n = nrow(sub),
          n_cases = sum(resp == 1),
          n_controls = sum(resp == 0)
        ) %>%
        select(Superpopulation, n, n_cases, n_controls, Score,
               auc, auc_lower, auc_upper, nagelkerke_R2)
    })
  
  return(results)
}

compute_or_two_scores <- function(
  data,
  score1_col,
  score2_col,
  case_col,
  ancestry_col,
  covariates = c("Age"),
  pops = NULL  # optional vector of ancestry groups to keep
) {
  
  # Subset populations if requested
  if (!is.null(pops)) {
    data <- data %>% filter(.data[[ancestry_col]] %in% pops)
  }
  
  or_table <- data %>%
    group_split(across(all_of(ancestry_col))) %>%
    map_dfr(function(sub) {
      ancestry <- unique(sub[[ancestry_col]])
      
      n_cases <- sum(sub[[case_col]] == 1, na.rm = TRUE)
      n_controls <- sum(sub[[case_col]] == 0, na.rm = TRUE)
      
      # Only fit model if both cases and controls exist
      if (length(unique(sub[[case_col]])) > 1) {
        
        # Build formula string dynamically: include covariates if provided
        formula_str <- paste(case_col, "~", score1_col, "+", score2_col)
        if (!is.null(covariates) && length(covariates) > 0) {
          formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
        }
        
        # Fit logistic regression model
        model <- glm(as.formula(formula_str), data = sub, family = binomial(link = "logit"))
        tidy_model <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
        
        # Extract ORs and p-values for the two scores
        scores_to_extract <- c(score1_col, score2_col)
        tidy_model %>%
          filter(term %in% scores_to_extract) %>%
          transmute(
            Superpopulation = ancestry,
            n = nrow(sub),
            n_cases = n_cases,
            n_controls = n_controls,
            Score = term,
            OR = round(estimate, 2),
            OR_lower = round(conf.low, 4),
            OR_upper = round(conf.high, 4),
            p_value = format.pval(p.value,
                                  digits = 4,
                                  scientific = TRUE,
                                  eps = 1E-300)
          )
        
      } else {
        # If only one outcome class
        tibble(
          Superpopulation = ancestry,
          n = nrow(sub),
          n_cases = n_cases,
          n_controls = n_controls,
          Score = c(score1_col, score2_col),
          OR = NA_real_,
          OR_lower = NA_real_,
          OR_upper = NA_real_,
          p_value = NA_character_
        )
      }
    })
  
  return(or_table)
}

compute_joint_or <- function(
  data,
  case_col,
  ancestry_col = NULL,
  CV_col = "CV_PGS_standardized",
  RV_col = "RV_PGS_standardized",
  CV_breaks = c(0, 0.4, 0.6, 1),
  RV_breaks = c(0, 0.1, 0.2, 0.8, 0.9, 1),
  CV_labels = NULL,
  RV_labels = NULL,
  ref_CV = 2,
  ref_RV = 3,
  covariates = NULL,
  pops = "EUR"
) {
  # Ensure unique breaks
  CV_breaks <- quantile(data[[CV_col]], probs = CV_breaks, na.rm = TRUE) %>% unique()
  RV_breaks <- quantile(data[[RV_col]], probs = RV_breaks, na.rm = TRUE) %>% unique()

  # Filter population
  data <- data %>%
    filter(Superpopulation %in% pops) %>%
    mutate(
      CV_group = cut(.data[[CV_col]], breaks = CV_breaks, include.lowest = TRUE, labels = FALSE),
      RV_group = cut(.data[[RV_col]], breaks = RV_breaks, include.lowest = TRUE, labels = FALSE)
    )

  nCV <- length(unique(data$CV_group))
  nRV <- length(unique(data$RV_group))

  # Compute joint order
  data <- data %>%
    mutate(JointOrder = (CV_group - 1) * nRV + RV_group)

  ref_joint <- (ref_CV - 1) * nRV + ref_RV

  # Labels
  if (is.null(CV_labels)) CV_labels <- paste0("CV", 1:nCV)
  if (is.null(RV_labels)) RV_labels <- paste0("RV", 1:nRV)

  data <- data %>%
    mutate(
      JointLabel = paste0(CV_labels[CV_group], ", ", RV_labels[RV_group])
    )

  # Factor levels: logical order by CV_group then RV_group
  joint_levels <- data %>%
    distinct(JointOrder, JointLabel, CV_group, RV_group) %>%
    arrange(CV_group, RV_group) %>%
    pull(JointLabel)

  data$JointLabel <- factor(data$JointLabel, levels = joint_levels)

  # Split by ancestry
  if (!is.null(ancestry_col)) {
    split_data <- data %>% group_split(.data[[ancestry_col]])
  } else {
    split_data <- list(data)
  }

  # Compute OR
  compute_one <- function(sub) {
    ancestry <- if (!is.null(ancestry_col)) unique(sub[[ancestry_col]]) else "All"

    counts_tbl <- sub %>%
      group_by(JointLabel) %>%
      summarise(
        n = n(),
        n_cases = sum(.data[[case_col]] == 1, na.rm = TRUE),
        n_controls = sum(.data[[case_col]] == 0, na.rm = TRUE),
        .groups = "drop"
      ) %>%
  mutate(
    n_cases = case_when(
      n_cases <= 20 ~ "≤ 20",
      TRUE ~ as.character(n_cases)
    )
  )

    # Handle no variation
    if (length(unique(sub[[case_col]])) < 2) {
      return(counts_tbl %>% mutate(OR=NA, OR_lower=NA, OR_upper=NA, p_value=NA, Superpopulation=ancestry))
    }

    # Logistic regression
    sub$JointOrder <- factor(sub$JointOrder, levels = sort(unique(sub$JointOrder)))
    sub$JointOrder <- relevel(sub$JointOrder, ref = as.character(ref_joint))

    formula_str <- paste(case_col, "~ JointOrder")
    if (!is.null(covariates)) formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))

    model <- glm(as.formula(formula_str), data = sub, family = binomial())
    tidy_m <- tidy(model, conf.int = TRUE, exponentiate = TRUE)

    df <- counts_tbl %>%
      mutate(OR = 1, OR_lower = 1, OR_upper = 1, p_value = NA_character_, Superpopulation = ancestry)

    for (i in seq_len(nrow(tidy_m))) {
      term <- tidy_m$term[i]
      if (term != "(Intercept)") {
        level_num <- as.numeric(gsub("JointOrder", "", term))
        joint_label <- sub$JointLabel[sub$JointOrder == level_num][1]
        df$OR[df$JointLabel == joint_label] <- round(tidy_m$estimate[i], 3)
        df$OR_lower[df$JointLabel == joint_label] <- round(tidy_m$conf.low[i], 3)
        df$OR_upper[df$JointLabel == joint_label] <- round(tidy_m$conf.high[i], 3)
        df$p_value[df$JointLabel == joint_label] <- format.pval(tidy_m$p.value[i], digits=3, eps=1e-300)
      }
    }

    df
  }

  map_dfr(split_data, compute_one)
}

plot_joint_or <- function(or_table) {

  ggplot(or_table, aes(x = x_lab, y = OR, ymin = OR_lower, ymax = OR_upper,
                       color = color_lab)) +
    geom_pointrange(position = position_dodge(width = 0.6), size = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
    labs(x = "CV PGS quantile", y = "OR", color = "RV PGS quantile") +
    scale_y_log10(breaks = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 10), minor_breaks = NULL) +
    theme_classic(base_size = 10) +
    scale_color_manual(
    values = c("darkgreen","black", "#CC3399", "#FF0000"),) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 14)
    ) + facet_wrap(~Superpopulation)
}


base_aes <- aes(color = Superpopulation, fill = Superpopulation)

# Compute shared axis limits
xlim <- range(breast_raw_data$SCORE1_AVG, na.rm = TRUE) * c(1.01,1.01)
ylim <- range(breast_raw_data$raw_RV_PGS, na.rm = TRUE) * c(1.01,1.01)

cor_val <- cor(
  breast_raw_data$SCORE1_AVG,
  breast_raw_data$raw_RV_PGS,
  use = "complete.obs",
  method = "pearson"
)

cor_label <- sprintf("Pearson r = %.3f", cor_val)

# Scatterplot (two scores)
p_scatter <- ggplot(breast_raw_data, aes(x = SCORE1_AVG,
                                  y = raw_RV_PGS,
                                  color = Superpopulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  annotate(
    "text",
    x = xlim[1] + 0.02 * diff(xlim),
    y = ylim[2] - 0.04 * diff(ylim),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "black"
  ) +

  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Raw CV PGS", y = "Raw RV PGS") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(breast_raw_data, aes(x = SCORE1_AVG, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(breast_raw_data, aes(x = raw_RV_PGS, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_breast_raw_PGS_values.jpg",plot,dpi = 800, height = 10, width = 10)

# CV PGS

original_names = names(breast_raw_data)

breast_standardized_data <- standardize_pgs_pca(
  breast_raw_data,
  pgs_col = "SCORE1_AVG",
  pca_cols = paste0("PC", 1:30),
  ref_filter = "breast_cancer == 0"
)


breast_standardized_data = breast_standardized_data %>% 
  select(all_of(original_names), Z_norm2)

colnames(breast_standardized_data) = c(original_names,"CV_PGS_standardized")

# RV PGS

## Look for quantile gap (largest relative jump between adjacent quantiles)
x <- breast_standardized_data$raw_RV_PGS

qs <- quantile(x, probs = seq(0.1, 0.5, by = 0.01))
ratios <- qs[-1] / qs[-length(qs)]

i <- which.max(ratios)

## Select threshold
threshold <- qs[i]

## Apply threshold to separate most RV carriers and non-carriers
breast_standardized_data <- breast_standardized_data %>%
  mutate(
    RV_group = ifelse(raw_RV_PGS <= threshold, "likely_non_carrier", "likely_carrier")
  )

## Apply square root transformation to RV PGS
breast_standardized_data <- breast_standardized_data %>%
  mutate(
    RV_PGS_standardized = sqrt(raw_RV_PGS)
  )

threshold

breast_standardized_data = breast_standardized_data %>%
    mutate(case = ifelse(breast_cancer == 1, "Case", "Control") %>% as.factor())

case_colors <- c(
  "Case"    = "red",
  "Control" = "black"
)

cor_val <- cor(
  breast_standardized_data$SCORE1_AVG,
  breast_standardized_data$raw_RV_PGS,
  use = "complete.obs",
  method = "pearson"
)

cor_label <- sprintf("Pearson r = %.3f", cor_val)


base_aes <- aes(color = case, fill = case)

# Compute shared axis limits
xlim <- range(breast_standardized_data$CV_PGS_standardized, na.rm = TRUE) * c(1.1,1.1)
ylim <- range(breast_standardized_data$RV_PGS_standardized, na.rm = TRUE) * c(0.99,1.1)
# Scatterplot (two scores)
p_scatter <- ggplot(
  breast_standardized_data,
  aes(x = CV_PGS_standardized, y = RV_PGS_standardized)
) +
  geom_point(
    data = subset(breast_standardized_data, case == "Control"),
    aes(color = case),
    alpha = 0.4, size = 1.2
  ) +
  annotate(
    "text",
    x = xlim[1] + 0.02 * diff(xlim),
    y = ylim[2] - 0.04 * diff(ylim),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "black"
  ) +

  geom_point(
    data = subset(breast_standardized_data, case == "Case"),
    aes(color = case),
    alpha = 0.7, size = 1.2
  ) +
  scale_color_manual(values = case_colors, name = "Breast cancer") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Standardized CV PGS", y = "sqrt(RV PGS)") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(
  breast_standardized_data,
  aes(x = CV_PGS_standardized, fill = case, color = case)
) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = case_colors, name = "Breast cancer") +
  scale_color_manual(values = case_colors, name = "Breast cancer") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(
  breast_standardized_data,
  aes(x = RV_PGS_standardized, fill = case, color = case)
) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = case_colors, name = "Breast cancer") +
  scale_color_manual(values = case_colors, name = "Breast cancer") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_breast_standardized_PGS_values_case_status.jpg",plot,dpi = 800, height = 10, width = 10)

base_aes <- aes(color = Superpopulation, fill = Superpopulation)

# Compute shared axis limits
xlim <- range(breast_standardized_data$CV_PGS_standardized, na.rm = TRUE) * c(1.1,1.1)
ylim <- range(breast_standardized_data$RV_PGS_standardized, na.rm = TRUE) * c(0.99,1.1)
# Scatterplot (two scores)
p_scatter <- ggplot(breast_standardized_data, aes(x = CV_PGS_standardized,
                                  y = RV_PGS_standardized,
                                  color = Superpopulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Standardized CV PGS", y = "sqrt(RV PGS)") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(breast_standardized_data, aes(x = CV_PGS_standardized, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(breast_standardized_data, aes(x = RV_PGS_standardized, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_breast_standardized_PGS_values.jpg",plot,dpi = 800, height = 10, width = 10)

summary_table <- evaluate_auc_r2_two_scores(
  data       = breast_standardized_data,
  score1_col = "CV_PGS_standardized",
  score2_col = "RV_PGS_standardized",
  case_col   = "breast_cancer",
  ancestry_col = "Superpopulation",
  pops = c("AFR","AMR","EUR")
)

  
summary_table

fwrite(summary_table, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_breast_internal_AUC.txt")

breast_standardized_data$sqrt_scaled_RV_PGS = 10*breast_standardized_data$RV_PGS_standardized

or_results <- compute_or_two_scores(
  data = breast_standardized_data,
  score1_col = "CV_PGS_standardized",
  score2_col = "sqrt_scaled_RV_PGS",
  case_col = "breast_cancer",
  ancestry_col = "Superpopulation",
  covariates = c("Age"),
  pops = c("AFR","AMR","EUR")
)

or_results
fwrite(or_results, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_breast_internal_OR_by_PGS_SD.txt")

threshold_name_percent = names(threshold)
threshold_decimal <- gsub("%", "", threshold_name_percent) %>% as.numeric() / 100

or_table <- compute_joint_or(
  data = breast_standardized_data,
  ancestry_col = "Superpopulation",
  case_col = "breast_cancer",
  CV_col = "CV_PGS_standardized",
  RV_col = "RV_PGS_standardized",
  CV_breaks = c(0, 0.33, 0.66, 1),
  RV_breaks = c(0, threshold_decimal, 0.75, 0.90, 1),
  CV_labels = c("[0–33%)", "[33–66%)", "[66–100%]"),
  RV_labels = c(
    paste0("[0–", threshold_name_percent, ")"),
    paste0("[", threshold_name_percent, "–75%)"),
    "[75–90%)",
    "[90–100%]"
  ),
  ref_CV = 2,
  ref_RV = 1,
  covariates = c("Age"),
  pops = c("AFR", "AMR", "EUR")
)

or_table

fwrite(or_table, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_breast_internal_OR_by_quantile_table.txt")

or_table = or_table %>%
  mutate(color_lab = sub("^[^,]*,\\s*", "", JointLabel),
         x_lab = sub(",.*$", "", JointLabel))


plot <- plot_joint_or(or_table)

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_breast_internal_OR_by_quantile_plot.jpg",plot,dpi=800,height=10,width=10)

base_aes <- aes(color = Superpopulation, fill = Superpopulation)

# Compute shared axis limits
xlim <- range(colorectal_raw_data$SCORE1_AVG, na.rm = TRUE) * c(1.01,1.01)
ylim <- range(colorectal_raw_data$raw_RV_PGS, na.rm = TRUE) * c(0.98,1.01)

cor_val <- cor(
  colorectal_raw_data$SCORE1_AVG,
  colorectal_raw_data$raw_RV_PGS,
  use = "complete.obs",
  method = "pearson"
)

cor_label <- sprintf("Pearson r = %.3f", cor_val)


# Scatterplot (two scores)
p_scatter <- ggplot(colorectal_raw_data, aes(x = SCORE1_AVG,
                                  y = raw_RV_PGS,
                                  color = Superpopulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  annotate(
    "text",
    x = xlim[1] + 0.02 * diff(xlim),
    y = ylim[2] - 0.04 * diff(ylim),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "black"
  ) +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Raw CV PGS", y = "Raw RV PGS") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(colorectal_raw_data, aes(x = SCORE1_AVG, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(colorectal_raw_data, aes(x = raw_RV_PGS, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_colorectal_raw_PGS_values.jpg",plot,dpi = 800, height = 10, width = 10)

# CV PGS

original_names = names(colorectal_raw_data)

colorectal_standardized_data <- standardize_pgs_pca(
  colorectal_raw_data,
  pgs_col = "SCORE1_AVG",
  pca_cols = paste0("PC", 1:30),
  ref_filter = "colorectal_cancer == 0"
)


colorectal_standardized_data = colorectal_standardized_data %>% 
  select(all_of(original_names), Z_norm2)

colnames(colorectal_standardized_data) = c(original_names,"CV_PGS_standardized")

# RV PGS

## Look for quantile gap (largest relative jump between adjacent quantiles)
x <- colorectal_standardized_data$raw_RV_PGS

qs <- quantile(x, probs = seq(0.1, 0.5, by = 0.01))
ratios <- qs[-1] / qs[-length(qs)]

i <- which.max(ratios)

## Select threshold
threshold <- qs[i]

## Apply threshold to separate most RV carriers and non-carriers
colorectal_standardized_data <- colorectal_standardized_data %>%
  mutate(
    RV_group = ifelse(raw_RV_PGS <= threshold, "likely_non_carrier", "likely_carrier")
  )

## Apply square root transformation to RV PGS
colorectal_standardized_data <- colorectal_standardized_data %>%
  mutate(
    RV_PGS_standardized = sqrt(raw_RV_PGS)
  )

colorectal_standardized_data = colorectal_standardized_data %>%
    mutate(case = ifelse(colorectal_cancer == 1, "Case", "Control") %>% as.factor())

case_colors <- c(
  "Case"    = "red",
  "Control" = "black"
)

base_aes <- aes(color = case, fill = case)

cor_val <- cor(
  colorectal_standardized_data$SCORE1_AVG,
  colorectal_standardized_data$raw_RV_PGS,
  use = "complete.obs",
  method = "pearson"
)

cor_label <- sprintf("Pearson r = %.3f", cor_val)


# Compute shared axis limits
xlim <- range(colorectal_standardized_data$CV_PGS_standardized, na.rm = TRUE) * c(1.1,1.1)
ylim <- range(colorectal_standardized_data$RV_PGS_standardized, na.rm = TRUE) * c(0.99,1.1)
# Scatterplot (two scores)
p_scatter <- ggplot(
  colorectal_standardized_data,
  aes(x = CV_PGS_standardized, y = RV_PGS_standardized)
) +
  geom_point(
    data = subset(colorectal_standardized_data, case == "Control"),
    aes(color = case),
    alpha = 0.4, size = 1.2
  ) +
  annotate(
    "text",
    x = xlim[1] + 0.02 * diff(xlim),
    y = ylim[2] - 0.04 * diff(ylim),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "black"
  ) +
  geom_point(
    data = subset(colorectal_standardized_data, case == "Case"),
    aes(color = case),
    alpha = 0.7, size = 1.2
  ) +
  scale_color_manual(values = case_colors, name = "colorectal cancer") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Standardized CV PGS", y = "sqrt(RV PGS)") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(
  colorectal_standardized_data,
  aes(x = CV_PGS_standardized, fill = case, color = case)
) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = case_colors, name = "colorectal cancer") +
  scale_color_manual(values = case_colors, name = "colorectal cancer") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(
  colorectal_standardized_data,
  aes(x = RV_PGS_standardized, fill = case, color = case)
) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = case_colors, name = "colorectal cancer") +
  scale_color_manual(values = case_colors, name = "colorectal cancer") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_colorectal_standardized_PGS_values_case_status.jpg",plot,dpi = 800, height = 10, width = 10)

base_aes <- aes(color = Superpopulation, fill = Superpopulation)

# Compute shared axis limits
xlim <- range(colorectal_standardized_data$CV_PGS_standardized, na.rm = TRUE) * c(1.1,1.1)
ylim <- range(colorectal_standardized_data$RV_PGS_standardized, na.rm = TRUE) * c(0.99,1.1)
# Scatterplot (two scores)
p_scatter <- ggplot(colorectal_standardized_data, aes(x = CV_PGS_standardized,
                                  y = RV_PGS_standardized,
                                  color = Superpopulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Standardized CV PGS", y = "sqrt(RV PGS)") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(colorectal_standardized_data, aes(x = CV_PGS_standardized, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(colorectal_standardized_data, aes(x = RV_PGS_standardized, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_colorectal_standardized_PGS_values.jpg",plot,dpi = 800, height = 10, width = 10)

summary_table <- evaluate_auc_r2_two_scores(
  data       = colorectal_standardized_data,
  score1_col = "CV_PGS_standardized",
  score2_col = "RV_PGS_standardized",
  case_col   = "colorectal_cancer",
  ancestry_col = "Superpopulation",
  pops = c("AFR","AMR","EUR")
)

  
summary_table

fwrite(summary_table, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_colorectal_internal_AUC.txt")


colorectal_standardized_data$sqrt_scaled_RV_PGS = 10*colorectal_standardized_data$RV_PGS_standardized

or_results <- compute_or_two_scores(
  data = colorectal_standardized_data,
  score1_col = "CV_PGS_standardized",
  score2_col = "sqrt_scaled_RV_PGS",
  case_col = "colorectal_cancer",
  ancestry_col = "Superpopulation",
  covariates = c("Age","Female"),
  pops = c("AFR","AMR","EUR")
)

or_results
fwrite(or_results, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_colorectal_internal_OR_by_PGS_SD.txt")

threshold_name_percent = names(threshold)
threshold_decimal <- gsub("%", "", threshold_name_percent) %>% as.numeric() / 100

or_table <- compute_joint_or(
  data = colorectal_standardized_data,
  ancestry_col = "Superpopulation",
  case_col = "colorectal_cancer",
  CV_col = "CV_PGS_standardized",
  RV_col = "RV_PGS_standardized",
  CV_breaks = c(0, 0.33, 0.66, 1),
  RV_breaks = c(0, threshold_decimal, 0.75, 0.90, 1),
  CV_labels = c("[0–33%)", "[33–66%)", "[66–100%]"),
  RV_labels = c(
    paste0("[0–", threshold_name_percent, ")"),
    paste0("[", threshold_name_percent, "–75%)"),
    "[75–90%)",
    "[90–100%]"
  ),
  ref_CV = 2,
  ref_RV = 1,
  covariates = c("Age","Female"),
  pops = c("AFR", "AMR", "EUR")
)

or_table

fwrite(or_table, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_colorectal_internal_OR_by_quantile_table.txt")

or_table = or_table %>%
  mutate(color_lab = sub("^[^,]*,\\s*", "", JointLabel),
         x_lab = sub(",.*$", "", JointLabel))


plot <- plot_joint_or(or_table)

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_colorectal_internal_OR_by_quantile_plot.jpg",plot,dpi=800,height=10,width=10)

base_aes <- aes(color = Superpopulation, fill = Superpopulation)

# Compute shared axis limits
xlim <- range(prostate_raw_data$SCORE1_AVG, na.rm = TRUE) * c(0.99,1.01)
ylim <- range(prostate_raw_data$raw_RV_PGS, na.rm = TRUE) * c(1.01,1.01)

cor_val <- cor(
  prostate_raw_data$SCORE1_AVG,
  prostate_raw_data$raw_RV_PGS,
  use = "complete.obs",
  method = "pearson"
)

cor_label <- sprintf("Pearson r = %.3f", cor_val)

# Scatterplot (two scores)
p_scatter <- ggplot(prostate_raw_data, aes(x = SCORE1_AVG,
                                  y = raw_RV_PGS,
                                  color = Superpopulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  annotate(
    "text",
    x = xlim[1] + 0.02 * diff(xlim),
    y = ylim[2] - 0.04 * diff(ylim),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "black"
  ) +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Raw CV PGS", y = "Raw RV PGS") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(prostate_raw_data, aes(x = SCORE1_AVG, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(prostate_raw_data, aes(x = raw_RV_PGS, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_prostate_raw_PGS_values.jpg",plot,dpi = 800, height = 10, width = 10)

# CV PGS

original_names = names(prostate_raw_data)

prostate_standardized_data <- standardize_pgs_pca(
  prostate_raw_data,
  pgs_col = "SCORE1_AVG",
  pca_cols = paste0("PC", 1:30),
  ref_filter = "prostate_cancer == 0"
)


prostate_standardized_data = prostate_standardized_data %>% 
  select(all_of(original_names), Z_norm2)

colnames(prostate_standardized_data) = c(original_names,"CV_PGS_standardized")

# RV PGS

## Look for quantile gap (largest relative jump between adjacent quantiles)
x <- prostate_standardized_data$raw_RV_PGS

qs <- quantile(x, probs = seq(0.1, 0.5, by = 0.01))
ratios <- qs[-1] / qs[-length(qs)]

i <- which.max(ratios)

## Select threshold
threshold <- qs[i]

## Apply threshold to separate most RV carriers and non-carriers
prostate_standardized_data <- prostate_standardized_data %>%
  mutate(
    RV_group = ifelse(raw_RV_PGS <= threshold, "likely_non_carrier", "likely_carrier")
  )

## Apply square root transformation to RV PGS
prostate_standardized_data <- prostate_standardized_data %>%
  mutate(
    RV_PGS_standardized = sqrt(raw_RV_PGS)
  )


prostate_standardized_data = prostate_standardized_data %>%
    mutate(case = ifelse(prostate_cancer == 1, "Case", "Control") %>% as.factor())

case_colors <- c(
  "Case"    = "red",
  "Control" = "black"
)

cor_val <- cor(
  prostate_standardized_data$SCORE1_AVG,
  prostate_standardized_data$raw_RV_PGS,
  use = "complete.obs",
  method = "pearson"
)

cor_label <- sprintf("Pearson r = %.3f", cor_val)

base_aes <- aes(color = case, fill = case)

# Compute shared axis limits
xlim <- range(prostate_standardized_data$CV_PGS_standardized, na.rm = TRUE) * c(1.1,1.1)
ylim <- range(prostate_standardized_data$RV_PGS_standardized, na.rm = TRUE) * c(0.99,1.1)
# Scatterplot (two scores)
p_scatter <- ggplot(
  prostate_standardized_data,
  aes(x = CV_PGS_standardized, y = RV_PGS_standardized)
) +
  geom_point(
    data = subset(prostate_standardized_data, case == "Control"),
    aes(color = case),
    alpha = 0.4, size = 1.2
  ) +
  annotate(
    "text",
    x = xlim[1] + 0.02 * diff(xlim),
    y = ylim[2] - 0.04 * diff(ylim),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    color = "black"
  ) +
  geom_point(
    data = subset(prostate_standardized_data, case == "Case"),
    aes(color = case),
    alpha = 0.7, size = 1.2
  ) +
  scale_color_manual(values = case_colors, name = "prostate cancer") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Standardized CV PGS", y = "sqrt(RV PGS)") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(
  prostate_standardized_data,
  aes(x = CV_PGS_standardized, fill = case, color = case)
) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = case_colors, name = "prostate cancer") +
  scale_color_manual(values = case_colors, name = "prostate cancer") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(
  prostate_standardized_data,
  aes(x = RV_PGS_standardized, fill = case, color = case)
) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_manual(values = case_colors, name = "prostate cancer") +
  scale_color_manual(values = case_colors, name = "prostate cancer") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_prostate_standardized_PGS_values_case_status.jpg",plot,dpi = 800, height = 10, width = 10)

base_aes <- aes(color = Superpopulation, fill = Superpopulation)

# Compute shared axis limits
xlim <- range(prostate_standardized_data$CV_PGS_standardized, na.rm = TRUE) * c(1.1,1.1)
ylim <- range(prostate_standardized_data$RV_PGS_standardized, na.rm = TRUE) * c(0.99,1.1)
# Scatterplot (two scores)
p_scatter <- ggplot(prostate_standardized_data, aes(x = CV_PGS_standardized,
                                  y = RV_PGS_standardized,
                                  color = Superpopulation)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  scale_y_continuous(limits = ylim, expand = c(0,0)) +
  labs(x = "Standardized CV PGS", y = "sqrt(RV PGS)") +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for X
p_xdens <- ggplot(prostate_standardized_data, aes(x = CV_PGS_standardized, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = xlim, expand = c(0,0)) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

# Density plot for Y (horizontal)
p_ydens <- ggplot(prostate_standardized_data, aes(x = RV_PGS_standardized, fill = Superpopulation, color = Superpopulation)) +
  geom_density(alpha = 0.4, adjust = 1.2) +
  scale_fill_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_color_brewer(palette = "Set2", name = "Inferred 1000G Superpopulation") +
  scale_x_continuous(limits = ylim, expand = c(0,0)) +
  coord_flip() +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0,0,0,0)
  )

layout <- "
AB
CD
"

# A = top density (x)
# B = empty space
# C = scatter
# D = right density

plot = (p_xdens + patchwork::plot_spacer() +
  p_scatter + p_ydens) +
  plot_layout(
    design = layout,
    heights = c(2, 4),
    widths = c(4, 2)
  )

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_internal_prostate_standardized_PGS_values.jpg",plot,dpi = 800, height = 10, width = 10)

summary_table <- evaluate_auc_r2_two_scores(
  data       = prostate_standardized_data,
  score1_col = "CV_PGS_standardized",
  score2_col = "RV_PGS_standardized",
  case_col   = "prostate_cancer",
  ancestry_col = "Superpopulation",
  pops = c("AFR","AMR","EUR")
)

  
summary_table

fwrite(summary_table, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_prostate_internal_AUC.txt")


prostate_standardized_data$sqrt_scaled_RV_PGS = 10*prostate_standardized_data$RV_PGS_standardized

or_results <- compute_or_two_scores(
  data = prostate_standardized_data,
  score1_col = "CV_PGS_standardized",
  score2_col = "sqrt_scaled_RV_PGS",
  case_col = "prostate_cancer",
  ancestry_col = "Superpopulation",
  covariates = c("Age"),
  pops = c("AFR","AMR","EUR")
)

or_results
fwrite(or_results, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_prostate_internal_OR_by_PGS_SD.txt")

threshold_name_percent = names(threshold)
threshold_decimal <- gsub("%", "", threshold_name_percent) %>% as.numeric() / 100

or_table <- compute_joint_or(
  data = prostate_standardized_data,
  ancestry_col = "Superpopulation",
  case_col = "prostate_cancer",
  CV_col = "CV_PGS_standardized",
  RV_col = "RV_PGS_standardized",
  CV_breaks = c(0, 0.33, 0.66, 1),
  RV_breaks = c(0, threshold_decimal, 0.75, 0.90, 1),
  CV_labels = c("[0–33%)", "[33–66%)", "[66–100%]"),
  RV_labels = c(
    paste0("[0–", threshold_name_percent, ")"),
    paste0("[", threshold_name_percent, "–75%)"),
    "[75–90%)",
    "[90–100%]"
  ),
  ref_CV = 2,
  ref_RV = 1,
  covariates = c("Age"),
  pops = c("AFR", "AMR", "EUR")
)

or_table

fwrite(or_table, "/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_prostate_internal_OR_by_quantile_table.txt")

or_table = or_table %>%
  mutate(color_lab = sub("^[^,]*,\\s*", "", JointLabel),
         x_lab = sub(",.*$", "", JointLabel))


plot <- plot_joint_or(or_table)

plot

ggplot2::ggsave("/home/jupyter/workspaces/v8geneticepidemiologyofcommoncancersinuspopulations/Aim_3/result_summaries/AoU_trained/AoU_trained_prostate_internal_OR_by_quantile_plot.jpg",plot,dpi=800,height=10,width=10)


