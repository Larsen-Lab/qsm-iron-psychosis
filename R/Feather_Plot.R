feather_plot <- function(df,
                         age.col   = "age",
                         id.col    = "bblid",
                         color.col = NA,
                         shape.col = NA,
                         color.values = NULL,
                         shape.values = NULL,
                         color.labels = waiver(),
                         shape.labels = waiver(),
                         x.breaks  = 5,
                         legend.position = c(0.85, 0.18)) {
  
  df <- df %>%
    mutate(id  = factor(.data[[id.col]]),
           age = .data[[age.col]]) %>%
    filter(!is.na(age))
  
  sorted.df <- df %>%
    group_by(id) %>%
    summarise(m = min(age, na.rm = TRUE), .groups = "drop") %>%
    arrange(m) %>%
    mutate(row = row_number())
  
  newdf <- df %>% left_join(sorted.df %>% select(id, row), by = "id")
  
  if (!is.na(color.col)) newdf$color.var <- factor(newdf[[color.col]])
  if (!is.na(shape.col)) newdf$shape.var <- factor(newdf[[shape.col]])
  
  base_aes <- aes(x = row, y = age, group = id)
  if (!is.na(color.col)) base_aes$colour <- quo(color.var)
  if (!is.na(shape.col)) base_aes$shape  <- quo(shape.var)
  
  if (is.null(color.values) && !is.na(color.col)) {
    okabe <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
               "#E69F00", "#56B4E9", "#F0E442", "#000000")
    lvls  <- levels(newdf$color.var)
    color.values <- setNames(okabe[seq_along(lvls)], lvls)
  }
  if (is.null(shape.values) && !is.na(shape.col)) {
    shp  <- c(16, 17, 15, 18, 8, 4, 3, 7)
    lvls <- levels(newdf$shape.var)
    shape.values <- setNames(shp[seq_along(lvls)], lvls)
  }
  
  n_part <- nrow(sorted.df)
  brks   <- unique(round(seq(1, n_part, length.out = x.breaks)))
  
  p <- ggplot(newdf, base_aes) +
    geom_line(alpha = 0.5, linewidth = 0.4) +
    geom_point(size = 1.2, alpha = 0.8, stroke = 0.3) +
    scale_x_continuous(breaks = brks, limits = c(1, n_part)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    coord_flip(clip = "off") +
    labs(x = "Participants",
         y = "Age (years)",
         color = NULL, shape = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.line   = element_line(linewidth = 0.35, color = "grey20"),
      axis.ticks  = element_line(linewidth = 0.35, color = "grey20"),
      legend.position   = legend.position,
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.key        = element_blank(),
      legend.key.size   = unit(0.4, "cm"),
      legend.spacing.y  = unit(2, "pt"),
      legend.margin     = margin(2, 4, 2, 4)
    )
  
  if (!is.na(color.col)) {
    p <- p + scale_color_manual(values = color.values, labels = color.labels) +
      guides(color = guide_legend(override.aes = list(linewidth = 1.5, size = 3)))
  }
  if (!is.na(shape.col)) {
    p <- p + scale_shape_manual(values = shape.values, labels = shape.labels) +
      guides(shape = guide_legend(override.aes = list(size = 3)))
  }
  
  p
}