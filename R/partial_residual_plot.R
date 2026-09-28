# Helper function for confidence intervals
compute_combined_se_gam <- function(modobj, terms_to_combine, newdata = NULL) {
  
  # A. Handle GAMM vs GAM objects
  if ("gamm" %in% class(modobj)) {
    modobj <- modobj$gam
  }
  
  # B. Get linear predictor matrix (Xp)
  # IMPORTANT: Use newdata if provided (for plotting curves), otherwise use model data
  if (!is.null(newdata)) {
    Xp <- predict(modobj, newdata = newdata, type = "lpmatrix")
  } else {
    Xp <- predict(modobj, type = "lpmatrix")
  }
  
  # C. Get variance-covariance matrix
  Vp <- vcov(modobj, unconditional = TRUE)
  
  # D. Identify which columns correspond to the terms
  smooth_labels <- sapply(modobj$smooth, function(x) x$label)
  param_names   <- names(coef(modobj)) # Fixed effects
  
  # Collect column indices for all terms
  all_cols <- c()
  
  for (tm in terms_to_combine) {
    
    # 1. Search Smooth Terms
    # Use fixed() to avoid regex errors with parenthesis
    smooth_matches <- which(str_detect(smooth_labels, fixed(tm)))
    
    if (length(smooth_matches) > 0) {
      for (i in smooth_matches) {
        smooth_obj <- modobj$smooth[[i]]
        all_cols <- c(all_cols, smooth_obj$first.para:smooth_obj$last.para)
      }
    }
    
    # 2. Search Parametric Terms (Fixed Effects)
    param_matches <- which(str_detect(param_names, fixed(tm)))
    
    if (length(param_matches) > 0) {
      all_cols <- c(all_cols, param_matches)
    }
  }
  
  all_cols <- unique(all_cols)
  
  # E. Safety Check
  if (length(all_cols) == 0) {
    warning(paste("Could not find terms:", paste(terms_to_combine, collapse=", "), 
                  "- Returning NULL"))
    return(NULL)
  }
  
  # F. Extract relevant sub-matrices
  Xp_sub <- Xp[, all_cols, drop = FALSE]
  Vp_sub <- Vp[all_cols, all_cols, drop = FALSE]
  
  # G. Compute SE: sqrt(diag(X * V * X^T))
  se_combined <- sqrt(rowSums((Xp_sub %*% Vp_sub) * Xp_sub))
  
  return(se_combined)
}

partial_resid_plot <- function(modobj, term, int.term = NULL, group.term = NULL, add.intercept = FALSE,
                               xlabel = NULL, ylabel = NULL, plot.title = NULL, show.data = FALSE,
                               newdata = NULL, filter.val = NULL,
                               focus_context = FALSE, diverge_age = NULL) {
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(mgcv)
  
  terms <- c(term, int.term)
  
  # Identify model type and extract data
  if ("gam" %in% class(modobj)) {
    model_type <- "gam"
    df <- if (!is.null(newdata)) newdata else modobj$model
  } else if ("gamm" %in% class(modobj) | "list" %in% class(modobj)) {
    model_type <- "gamm"
    df <- if (!is.null(newdata)) newdata else modobj$gam$model  # Extract model from the 'gam' component
    modobj <- modobj$gam  # Set the modobj to the 'gam' component for plotting
  } else if ("lm" %in% class(modobj)) {
    model_type <- "lm"
    df <- if (!is.null(newdata)) newdata else model.frame(modobj)
  } else if ("lmerMod" %in% class(modobj)) {
    model_type <- "lmer"
    df <- if (!is.null(newdata)) newdata else modobj@frame  # Use @frame for lmerMod objects
  } else {
    stop("Unsupported model type. Use gam, gamm, lm, or lmer.")
  }
  
  # Extract intercept if present
  mod.intercept <- if (add.intercept) {
    if (!"(Intercept)" %in% names(fixef(modobj))) {
      warning("Intercept not found in model coefficients; setting to 0.")
      0
    } else {
      fixef(modobj)["(Intercept)"]
    }
  } else {
    0
  }
  
  # Handling gam or gamm models
  if (model_type == "gam" || model_type == "gamm") {
    pterms <- predict(modobj, type = "terms", se.fit = TRUE, unconditional = TRUE)

    if (!is.null(attr(pterms$fit, "constant"))) {
      pterms$fit <- sweep(pterms$fit, 2, attr(pterms$fit, "constant"), "+")
    }

    # Add RE smooth contributions back to partial residuals so that between-subject
    # scatter is preserved in the display (otherwise conditional residuals are tiny)
    re_contribution <- {
      is_re <- sapply(modobj$smooth, function(s) "random.effect" %in% class(s))
      if (any(is_re)) {
        re_labels <- sapply(modobj$smooth[is_re], function(s) s$label)
        re_cols   <- which(colnames(pterms$fit) %in% re_labels)
        if (length(re_cols) > 0) rowSums(pterms$fit[, re_cols, drop = FALSE]) else rep(0, nrow(pterms$fit))
      } else {
        rep(0, nrow(pterms$fit))
      }
    }
    
    # Handle interaction terms if specified
    if (!is.null(int.term)) {
      intcol <- df %>% dplyr::select(int.term)
      int_type <- ifelse(any(class(intcol[[1]]) == "numeric"), "numeric", "factor")
      if (int_type=="numeric") {
        
      }
    }
    
    # Adjust for intercept if needed
    if (add.intercept) {
      pterms.fit <- pterms$fit + mod.intercept
    } else {
      pterms.fit <- pterms$fit
    }
    
    pterms.sefit <- pterms$se.fit
    
    # Clean term names (e.g., removing "s(" and ")")
    colnames(pterms.fit) <- gsub(x = colnames(pterms.fit), pattern = "s\\(", replacement = "") %>%
      gsub(pattern = "\\)", replacement = "")
    colnames(pterms.sefit) <- gsub(x = colnames(pterms.sefit), pattern = "s\\(", replacement = "") %>%
      gsub(pattern = "\\)", replacement = "")
    
    if (!is.null(int.term) && length(terms) > 1) {
      #### NOTE ####
      # Using variance-covariance matrix to compute combined SEs
      # This accounts for correlation between terms (e.g., main effect + interaction)
      ##############
      
      se_combined <- compute_combined_se_gam(modobj, terms_to_combine = terms, newdata = df)
      
      if (!is.null(se_combined)) {
        # Successfully computed exact SEs
        pterms.df <- data.frame(pterms.fit) %>%
          dplyr::select(matches(terms)) %>%
          mutate(fit = rowSums(across(where(is.numeric)))) %>%
          mutate(se.fit = se_combined,
                 upr = fit + 1.96 * se_combined,
                 lwr = fit - 1.96 * se_combined) %>%
          mutate(partial.residuals = fit + resid(modobj) + re_contribution) %>%
          dplyr::select(fit, se.fit, lwr, upr, partial.residuals)
      } else {
        # Fall back to approximate method
        warning("Falling back to approximate SE calculation (assumes independence)")
        pterms.df <- data.frame(pterms.fit) %>%
          dplyr::select(matches(terms)) %>%
          mutate(fit = rowSums(across(where(is.numeric)))) %>%
          dplyr::select(fit) %>%
          cbind(data.frame(pterms.sefit) %>%
                  dplyr::select(matches(terms)) %>%
                  mutate(se.fit = sqrt(rowSums(across(where(is.numeric))^2)))) %>%
          mutate(upr = fit + 1.96 * se.fit,
                 lwr = fit - 1.96 * se.fit) %>%
          mutate(partial.residuals = fit + resid(modobj) + re_contribution) %>%
          dplyr::select(fit, se.fit, lwr, upr, partial.residuals)
      }
    } else {
      # Single term - use standard approach
      pterms.df <- data.frame(pterms.fit) %>%
        dplyr::select(matches(term)) %>%
        mutate(fit = rowSums(across(where(is.numeric)))) %>%
        dplyr::select(fit) %>%
        cbind(data.frame(pterms.sefit) %>%
                dplyr::select(matches(term)) %>%
                mutate(se.fit = rowSums(across(where(is.numeric))))) %>%
        mutate(upr = fit + 1.96 * se.fit,
               lwr = fit - 1.96 * se.fit) %>%
        mutate(partial.residuals = fit + resid(modobj) + re_contribution) %>%
        dplyr::select(fit, se.fit, lwr, upr, partial.residuals)
    }
    
    partial.residuals.df <- cbind(pterms.df, df %>%
                                    dplyr::select(matches(term)) %>%
                                    rename_all(.funs = function(x) paste0(x, ".raw")))
    
    # Extract residuals and fits
    pterms <- partial.residuals.df
  } 
  
  # For lm models, extract terms and adjust for intercept if needed, and calculate residuals
  else if (model_type == "lm") {
    terms_pred <- predict(modobj, type = "terms", se.fit = TRUE)
    pterms.fit <- as.data.frame(terms_pred$fit)
    pterms.sefit <- as.data.frame(terms_pred$se.fit)
    
    if (add.intercept) {
      pterms.fit <- pterms + mod.intercept
    }
    
    # Compute residuals
    residuals <- resid(modobj)
    pterms.df <- data.frame(pterms.fit) %>%
      dplyr::select(matches(terms)) %>%
      mutate(fit = rowSums(across(where(is.numeric)))) %>%
      dplyr::select(fit) %>%
      cbind(data.frame(pterms.sefit) %>%
              dplyr::select(matches(terms)) %>%
              mutate(se.fit = rowSums(across(where(is.numeric))))) %>%
      mutate(upr = fit + 1.96 * se.fit,
             lwr = fit - 1.96 * se.fit) %>%
      mutate(partial.residuals = fit + residuals) %>%
      dplyr::select(fit, se.fit, lwr, upr, partial.residuals)
    pterms <- pterms.df
    
  } 
  
  # For lmer models, calculate residuals and combine with fixed-effects predictions for partial residuals
  else if (model_type == "lmer") {
    # Calculate fixed-effect contributions
    fixed_effects <- fixef(modobj)
    design_matrix <- model.matrix(terms(modobj), df)
    pterms <- as.data.frame(design_matrix %*% fixed_effects)  # Ensure pterms is a data frame
    
    # Add intercept if requested
    pterms <- pterms + if (add.intercept) mod.intercept else 0
    names(pterms) <- "fit"
    
    # Standard errors are not computed for lmer fixed effects in this implementation
    se_terms <- data.frame(rep(NA, nrow(df)))
    
    # Compute residuals
    residuals <- residuals(modobj)
    
    # Calculate partial residuals: fitted values + residuals
    pterms$partial.residuals <- pterms$fit + residuals
  } 
  
  else {
    stop("Unsupported model type.")
  }
  
  # Prepare plot data
  plot_data <- df %>%
    mutate(partial_residuals = pterms$partial.residuals,
           fit = pterms$fit, 
           lwr = if (model_type == "lmer") NA else pterms$lwr,  # Set NA for lwr and upr if model is lmer
           upr = if (model_type == "lmer") NA else pterms$upr) %>%
    dplyr::select(all_of(c(term, int.term, group.term)), partial_residuals, fit, lwr, upr)
  

  # Filter if needed (only if int.term is specified, since filter.val should correspond to levels of int.term)

  if (!is.null(filter.val) && !is.null(int.term)) {
    # Check if the interaction term exists in data
    if (int.term %in% names(plot_data)) {
      plot_data <- plot_data %>% 
        filter(.data[[int.term]] %in% filter.val)
      
      # Optional: Drop unused factor levels so they don't appear in legend
      if (is.factor(plot_data[[int.term]])) {
        plot_data[[int.term]] <- droplevels(plot_data[[int.term]])
      }
    } else {
      warning("filter.val was provided but int.term was not found in the data. No filtering applied.")
    }
  }

  
  # === Focus + Context mode ===
  if (focus_context && !is.null(int.term)) {

    has_focus <- !is.null(diverge_age) && is.numeric(diverge_age)

    # Alpha settings: muted context vs vivid focus
    a_ctx <- list(point = 0.3, spaghetti = 0.25, ribbon = 0.2, fit = 0.33)
    a_fcs <- list(point = 0.3, spaghetti = 0.25, ribbon = 0.45, fit = 1)

    # Build deduplicated fit curve (one row per unique age × group)
    fit_curve <- plot_data %>%
      group_by(across(all_of(c(term, int.term)))) %>%
      summarise(fit = mean(fit), lwr = mean(lwr), upr = mean(upr), .groups = "drop") %>%
      arrange(across(all_of(c(int.term, term))))

    if (has_focus) {
      # Interpolate boundary point per group level for seamless connection
      grps <- unique(fit_curve[[int.term]])
      bnd_rows <- lapply(grps, function(g) {
        sub <- fit_curve[fit_curve[[int.term]] == g, ]
        ages <- sub[[term]]
        if (diverge_age > min(ages) & diverge_age < max(ages) &
            !(diverge_age %in% ages)) {
          r <- sub[1, , drop = FALSE]
          r[[term]] <- diverge_age
          r$fit <- approx(ages, sub$fit, xout = diverge_age)$y
          r$lwr <- approx(ages, sub$lwr, xout = diverge_age)$y
          r$upr <- approx(ages, sub$upr, xout = diverge_age)$y
          r
        }
      })
      bnd_df <- do.call(rbind, bnd_rows)
      if (!is.null(bnd_df)) {
        fit_curve <- rbind(fit_curve, bnd_df) %>%
          arrange(across(all_of(c(int.term, term))))
      }

      ctx_fit <- fit_curve[fit_curve[[term]] <= diverge_age, ]
      fcs_fit <- fit_curve[fit_curve[[term]] >= diverge_age, ]
    } else {
      ctx_fit <- fit_curve
      fcs_fit <- fit_curve[0, ]
    }

    # Build plot
    gg <- ggplot(plot_data, aes(x = .data[[term]], y = partial_residuals))

    # Data points and spaghetti (uniformly muted)
    if (show.data) {
      gg <- gg +
        geom_point(aes(color = .data[[int.term]], fill = .data[[int.term]]),
                   alpha = a_ctx$point, shape = 21)
      if (!is.null(group.term)) {
        gg <- gg +
          geom_line(aes(group = .data[[group.term]], color = .data[[int.term]]),
                    alpha = a_ctx$spaghetti, size = 1)
      }
    }

    # Context ribbon + fit line (muted)
    if (model_type != "lmer") {
      gg <- gg +
        geom_ribbon(data = ctx_fit,
                    aes(x = .data[[term]], y = fit, ymin = lwr, ymax = upr,
                        fill = .data[[int.term]]),
                    alpha = a_ctx$ribbon, color = NA, inherit.aes = FALSE)
    }
    gg <- gg +
      geom_line(data = ctx_fit,
                aes(x = .data[[term]], y = fit, color = .data[[int.term]]),
                size = 1.5, alpha = a_ctx$fit, inherit.aes = FALSE)

    # Focus ribbon + fit line (vivid, post-divergence)
    if (has_focus && nrow(fcs_fit) > 0) {
      if (model_type != "lmer") {
        gg <- gg +
          geom_ribbon(data = fcs_fit,
                      aes(x = .data[[term]], y = fit, ymin = lwr, ymax = upr,
                          fill = .data[[int.term]]),
                      alpha = a_fcs$ribbon, color = NA, inherit.aes = FALSE)
      }
      gg <- gg +
        geom_line(data = fcs_fit,
                  aes(x = .data[[term]], y = fit, color = .data[[int.term]]),
                  size = 1.5, alpha = a_fcs$fit, inherit.aes = FALSE)

      # Subtle divergence marker with label
      gg <- gg +
        geom_vline(xintercept = diverge_age, linetype = "dashed",
                   color = "gray15", linewidth = .85) +
        annotate("text", x = diverge_age, y = -Inf,
                 label = paste0("Divergence point: ", round(diverge_age, 1), " yr"),
                 vjust = -0.5, hjust = -0.05, size = 4.5,
                 color = "gray15")
    }

    gg <- gg +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1")

  } else {
    # === Original plotting (unchanged) ===
    gg <- ggplot(plot_data, aes(x = .data[[term]], y = partial_residuals))

    if (show.data) {
      if (!is.null(int.term)) {
        gg <- gg + geom_point(aes(color = .data[[int.term]], fill = .data[[int.term]]),
                              alpha = 0.5, shape = 21)
        if (!is.null(group.term)) {
          gg <- gg + geom_line(aes(group = .data[[group.term]], color = .data[[int.term]]),
                               alpha = 0.45, size = 1)
        }
      } else {
        gg <- gg + geom_point(alpha = 0.5, color = "gray50", shape = 21, fill = "gray40", size = 1.5)
        if (!is.null(group.term)) {
          gg <- gg + geom_line(aes(group = .data[[group.term]]),
                               alpha = 0.45, size = .8, color = "gray10")
        }
      }
    }

    if (model_type != "lmer") {
      if (!is.null(int.term)) {
        gg <- gg + geom_ribbon(aes(y = fit, ymin = lwr, ymax = upr, fill = .data[[int.term]]),
                               alpha = 0.45, color = NA)
      } else {
        gg <- gg + geom_ribbon(aes(y = fit, ymin = lwr, ymax = upr),
                               fill = "gray70", alpha = 0.45, color = NA)
      }
    }

    if (!is.null(int.term)) {
      gg <- gg + geom_line(aes(y = fit, color = .data[[int.term]]),
                           size = 1.5, alpha = 1) +
        scale_fill_brewer(palette = "Set1") +
        scale_color_brewer(palette = "Set1")
    } else {
      gg <- gg + geom_line(aes(y = fit), color = "black", size = 1.5, alpha = 1)
    }
  }

  # Customize labels
  if (!is.null(xlabel)) gg <- gg + xlab(xlabel)
  if (!is.null(ylabel)) gg <- gg + ylab(ylabel)
  if (!is.null(plot.title)) gg <- gg + ggtitle(plot.title)

  return(gg)
}
