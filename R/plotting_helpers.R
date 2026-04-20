plot_manhattan_bh <- function(p, variants = NULL, k = 8, alpha = 0.1, DIR = "temp", file_name = "temp") {
  n <- length(p)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'p'")
  }
  
  # compute BH adjusted p-values
  padj <- p.adjust(p, method = "BH")
  
  # find BH significance threshold in raw p-value scale
  # (largest p that passes BH)
  bh_thresh <- max(p[padj <= alpha], na.rm = TRUE)
  
  # if nothing significant:
  if (is.infinite(bh_thresh)) bh_thresh <- NA
  
  df <- data.frame(
    pos = seq_along(p),
    logp = -log10(p),
    variant = if (!is.null(variants)) variants else NA
  )
  
  # Find k smallest p-values
  top_idx <- order(p)[1:min(k, n)]
  df$highlight <- FALSE
  df$highlight[top_idx] <- TRUE
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = logp)) +
    geom_point(aes(color = highlight), size = 1.5) +
    scale_color_manual(values = c("black", "blue")) +
    expand_limits(y = max(df$logp, na.rm = TRUE) * 1.02) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "-log10(p)",
      title = paste0("Manhattan Plot (BH FDR=", alpha, ")")
    ) +
    geom_hline(yintercept = -log10(bh_thresh), color = "red", linetype = "dashed") +
    theme(legend.position = "none") 
  
  # Add labels for top k
  if (!is.null(variants)) {
    plt <- plt +
      geom_text(
        data = df[top_idx, ],
        aes(label = variant),
        vjust = 1,
        size = 3
      )
  }
  
  # Make sure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  # Save
  ggsave(filename = file.path(DIR, file_name), plot = plt, width = 6, height = 4, units = "in",
         dpi = 300)
  
}

plot_manhattan_threshold <- function(p, variants = NULL, k = 8, threshold = 5e-08, DIR = "temp", file_name = "temp") {
  n <- length(p)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'p'")
  }
  
  df <- data.frame(
    pos = seq_along(p),
    logp = -log10(p),
    variant = if (!is.null(variants)) variants else NA
  )
  
  # Find k smallest p-values
  top_idx <- order(p)[1:min(k, n)]
  df$highlight <- FALSE
  df$highlight[top_idx] <- TRUE
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = logp)) +
    geom_point(aes(color = highlight), size = 1.5) +
    scale_color_manual(values = c("black", "blue")) +
    expand_limits(y = max(df$logp, na.rm = TRUE) * 1.02) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "-log10(p)",
      title = paste0("Manhattan Plot")
    ) + 
    geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed") +
    theme(legend.position = "none") 
  
  # Add labels for top k
  if (!is.null(variants)) {
    plt <- plt +
      geom_text(
        data = df[top_idx, ],
        aes(label = variant),
        vjust = 1,
        size = 3
      )
  }
  
  # Make sure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  # Save
  ggsave(filename = file.path(DIR, file_name), plot = plt, width = 6, height = 4, units = "in",
         dpi = 300)
  
}

plot_manhattan <- function(
    p, 
    variants = NULL, 
    k = 8, 
    alpha = 0.1,             # FDR for BH threshold
    threshold = 5e-08,        # user-specified p-value threshold (raw p)
    DIR = "temp", 
    file_name = "temp"
) {
  n <- length(p)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'p'")
  }
  
  # Filter out NA/Inf
  valid <- is.finite(p)
  p <- p[valid]
  if (!is.null(variants)) variants <- variants[valid]
  
  # BH threshold
  padj <- p.adjust(p, method = "BH")
  bh_thresh <- max(p[padj <= alpha], na.rm = TRUE)
  if (is.infinite(bh_thresh)) bh_thresh <- NA
  
  # Data frame
  df <- data.frame(
    pos = seq_along(p),
    logp = -log10(p),
    variant = if (!is.null(variants)) variants else NA
  )
  
  # Highlight top k smallest p-values
  top_idx <- order(p)[1:min(k, n)]
  df$highlight <- FALSE
  df$highlight[top_idx] <- TRUE
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = logp)) +
    geom_point(aes(color = highlight), size = 1.5) +
    scale_color_manual(values = c("black", "blue")) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "-log10(p)",
      title = paste0("Manhattan Plot (BH alpha = ", alpha, ", threshold = ", threshold, ")")
    ) + 
    theme(legend.position = "none") +
    expand_limits(y = max(df$logp, na.rm = TRUE) * 1.1)
  
  # Add horizontal lines for thresholds
  if (!is.na(bh_thresh)) {
    plt <- plt + geom_hline(yintercept = -log10(bh_thresh), color = "red", linetype = "dashed")
  }
  if (!is.null(threshold)) {
    plt <- plt + geom_hline(yintercept = -log10(threshold), color = "green", linetype = "dashed")
  }
  
  # Add labels for top k using ggrepel
  if (!is.null(variants)) {
    plt <- plt + 
      geom_text_repel(
        data = df[top_idx, ],
        aes(label = variant),
        nudge_y = 0.01,
        size = 3,
        max.overlaps = Inf
      )
  }
  
  # Make sure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  # Save
  ggsave(filename = file.path(DIR, file_name), plot = plt, width = 6, height = 4, units = "in",
         dpi = 300)
}

plot_manhattan_e <- function(
    e, 
    variants = NULL, 
    k = 8, 
    alpha = 0.1,             # FDR for eBH threshold
    pos_labels = NULL,
    neg_labels = NULL,
    DIR = "temp", 
    file_name = "temp"
) {
  n <- length(e)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'p'")
  }
  
  # Filter out NA/Inf
  valid <- is.finite(e)
  e <- e[valid]
  if (!is.null(variants)) variants <- variants[valid]
  
  # eBH threshold
  ebh_thresh <- eBH_threshold(e, alpha = alpha)
  if (is.infinite(ebh_thresh)) ebh_thresh <- NA
  
  # Data frame
  df <- data.frame(
    pos = seq_along(e),
    e = e,
    variant = if (!is.null(variants)) variants else NA
  )
  
  # Highlight top k smallest p-values
  top_idx <- order(-e)[1:min(k, n)]
  df$highlight <- FALSE
  df$highlight[top_idx] <- TRUE
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = e)) +
    geom_point(aes(color = highlight), size = 1.5) +
    scale_color_manual(values = c("black", "blue")) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "e",
      title = paste0("Manhattan Plot (eBH alpha = ", alpha, ")")
    ) + 
    theme(legend.position = "none") +
    expand_limits(y = max(df$e, na.rm = TRUE) * 1.1)
  
  # Add horizontal lines for thresholds
  if (!is.na(ebh_thresh)) {
    plt <- plt + geom_hline(yintercept = ebh_thresh, color = "red", linetype = "dashed")
  }
  
  # Add labels for top k using ggrepel
  if (!is.null(variants)) {
    plt <- plt + 
      geom_text_repel(
        data = df[top_idx, ],
        aes(label = variant),
        nudge_y = 0.01,
        size = 3,
        max.overlaps = Inf
      )
  }
  
  # Make sure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  # Save
  ggsave(filename = file.path(DIR, file_name), plot = plt, width = 6, height = 4, units = "in",
         dpi = 300)
}

eBH_threshold <- function(e_vals, alpha = 0.1) {
  e_vals <- as.numeric(e_vals)
  K <- length(e_vals)
  
  # sort from largest to smallest
  e_sorted <- sort(e_vals, decreasing = TRUE)
  
  # k = 1,...,K  (this is correct)
  k <- seq_len(K)
  
  # RHS = K/(alpha*k)
  rhs <- K / (alpha * k)
  
  valid <- which(e_sorted >= rhs)
  
  if (length(valid) == 0) {
    return(Inf)
  }
  
  k_star <- max(valid)
  return(e_sorted[k_star])
}

plot_manhattan_pip <- function(
    pip, 
    variants = NULL, 
    k = 8, 
    pos_labels = NULL,
    neg_labels = NULL,
    threshold = 0.9,
    DIR = "temp", 
    file_name = "temp"
) {
  n <- length(pip)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'pip'")
  }
  
  # Filter out NA/Inf
  valid <- is.finite(pip)
  pip <- pip[valid]
  if (!is.null(variants)) variants <- variants[valid]
  
  # Data frame
  df <- data.frame(
    pos = seq_along(pip),
    pip = pip,
    variant = if (!is.null(variants)) variants else NA,
    stringsAsFactors = FALSE
  )
  
  # ---- Assign category labels ----
  df$group <- "Unknown"
  
  if (!is.null(variants)) {
    if (!is.null(pos_labels) && length(pos_labels) > 0) {
      pos_match <- Reduce(`|`, lapply(pos_labels, function(s) grepl(s, df$variant, fixed = TRUE)))
      df$group[pos_match] <- "Confirmed +"
    }
    
    if (!is.null(neg_labels) && length(neg_labels) > 0) {
      neg_match <- Reduce(`|`, lapply(neg_labels, function(s) grepl(s, df$variant, fixed = TRUE)))
      df$group[neg_match] <- "Confirmed -"
    }
  }
  
  df$group <- factor(df$group, levels = c("Unknown", "Confirmed +", "Confirmed -"))
  
  # Top-k indices for labeling only
  top_idx <- order(-df$pip)[1:min(k, n)]
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = pip)) +
    geom_point(aes(color = group), size = 1.5) +
    scale_color_manual(
      values = c(
        "Unknown" = "black",
        "Confirmed +" = "green",
        "Confirmed -" = "red"
      ),
      name = "Variant class"
    ) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "pip",
      title = paste0("Manhattan Plot")
    ) +
    expand_limits(y = max(df$pip, na.rm = TRUE) * 1.1)
  
  # eBH threshold line
  if (!is.na(threshold)) {
    plt <- plt + 
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed")
  }
  
  # Labels for top k
  if (!is.null(variants)) {
    plt <- plt + 
      ggrepel::geom_text_repel(
        data = df[top_idx, ],
        aes(label = variant),
        size = 3,
        max.overlaps = Inf
      )
  }
  
  # Ensure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  ggsave(
    filename = file.path(DIR, file_name),
    plot = plt,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300
  )
}

plot_manhattan_e_test <- function(
    e, 
    variants = NULL, 
    k = 8, 
    alpha = 0.1,
    pos_labels = NULL,
    neg_labels = NULL,
    DIR = "temp", 
    file_name = "temp"
) {
  n <- length(e)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'e'")
  }
  
  # Filter out NA/Inf
  valid <- is.finite(e)
  e <- e[valid]
  if (!is.null(variants)) variants <- variants[valid]
  n <- length(e)
  
  # eBH threshold
  ebh_thresh <- eBH_threshold(e, alpha = alpha)
  if (is.infinite(ebh_thresh)) ebh_thresh <- NA
  
  # Data frame
  df <- data.frame(
    pos = seq_along(e),
    e = e,
    variant = if (!is.null(variants)) variants else NA,
    stringsAsFactors = FALSE
  )
  
  # ---- Assign category labels ----
  df$group <- "Unknown"
  
  if (!is.null(variants)) {
    if (!is.null(pos_labels) && length(pos_labels) > 0) {
      pos_match <- Reduce(`|`, lapply(pos_labels, function(s) grepl(s, df$variant, fixed = TRUE)))
      df$group[pos_match] <- "Confirmed +"
    }
    
    if (!is.null(neg_labels) && length(neg_labels) > 0) {
      neg_match <- Reduce(`|`, lapply(neg_labels, function(s) grepl(s, df$variant, fixed = TRUE)))
      df$group[neg_match] <- "Confirmed -"
    }
  }
  
  df$group <- factor(df$group, levels = c("Unknown", "Confirmed +", "Confirmed -"))
  
  # Top-k indices for labeling only
  top_idx <- order(-df$e)[1:min(k, n)]
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = e)) +
    geom_point(aes(color = group), size = 1.5) +
    scale_color_manual(
      values = c(
        "Unknown" = "black",
        "Confirmed +" = "green",
        "Confirmed -" = "red"
      ),
      name = "Variant class"
    ) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "e",
      title = paste0("Manhattan Plot (eBH alpha = ", alpha, ")")
    ) +
    expand_limits(y = max(df$e, na.rm = TRUE) * 1.1)
  
  # eBH threshold line
  if (!is.na(ebh_thresh)) {
    plt <- plt + 
      geom_hline(yintercept = ebh_thresh, color = "red", linetype = "dashed")
  }
  
  # Labels for top k
  if (!is.null(variants)) {
    plt <- plt + 
      ggrepel::geom_text_repel(
        data = df[top_idx, ],
        aes(label = variant),
        size = 3,
        max.overlaps = Inf
      )
  }
  
  # Ensure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  ggsave(
    filename = file.path(DIR, file_name),
    plot = plt,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300
  )
}

plot_manhattan_test <- function(
    p, 
    variants = NULL, 
    k = 8, 
    alpha = 0.1,
    threshold = 5e-08,
    pos_labels = NULL,
    neg_labels = NULL,
    DIR = "temp", 
    file_name = "temp"
) {
  n <- length(p)
  if (!is.null(variants) && length(variants) != n) {
    stop("Length of 'variants' must match length of 'p'")
  }
  
  # Filter out NA/Inf
  valid <- is.finite(p)
  p <- p[valid]
  if (!is.null(variants)) variants <- variants[valid]
  n <- length(p)
  
  # BH threshold
  padj <- p.adjust(p, method = "BH")
  bh_thresh <- suppressWarnings(max(p[padj <= alpha], na.rm = TRUE))
  if (is.infinite(bh_thresh)) bh_thresh <- NA
  
  # Data frame
  df <- data.frame(
    pos = seq_along(p),
    logp = -log10(p),
    variant = if (!is.null(variants)) variants else NA,
    stringsAsFactors = FALSE
  )
  
  # ---- Assign variant groups ----
  df$group <- "Unknown"
  
  if (!is.null(variants)) {
    if (!is.null(pos_labels) && length(pos_labels) > 0) {
      pos_match <- Reduce(`|`, lapply(pos_labels, function(s) grepl(s, df$variant, fixed = TRUE)))
      df$group[pos_match] <- "Confirmed +"
    }
    
    if (!is.null(neg_labels) && length(neg_labels) > 0) {
      neg_match <- Reduce(`|`, lapply(neg_labels, function(s) grepl(s, df$variant, fixed = TRUE)))
      df$group[neg_match] <- "Confirmed -"
    }
  }
  
  df$group <- factor(df$group, levels = c("Unknown", "Confirmed +", "Confirmed -"))
  
  # Top-k indices (for labeling only)
  top_idx <- order(p)[1:min(k, n)]
  
  # Plot
  plt <- ggplot(df, aes(x = pos, y = logp)) +
    geom_point(aes(color = group), size = 1.5) +
    scale_color_manual(
      values = c(
        "Unknown" = "black",
        "Confirmed +" = "green",
        "Confirmed -" = "red"
      ),
      name = "Variant class"
    ) +
    theme_minimal() +
    labs(
      x = "Index",
      y = "-log10(p)",
      title = paste0("Manhattan Plot (BH alpha = ", alpha, ", threshold = ", threshold, ")")
    ) +
    expand_limits(y = max(df$logp, na.rm = TRUE) * 1.1)
  
  # Threshold lines
  if (!is.na(bh_thresh)) {
    plt <- plt + 
      geom_hline(yintercept = -log10(bh_thresh), color = "blue", linetype = "dashed")
  }
  if (!is.null(threshold)) {
    plt <- plt + 
      geom_hline(yintercept = -log10(threshold), color = "lightblue", linetype = "dashed")
  }
  
  # Labels for top k
  if (!is.null(variants)) {
    plt <- plt + 
      ggrepel::geom_text_repel(
        data = df[top_idx, ],
        aes(label = variant),
        size = 3,
        max.overlaps = Inf
      )
  }
  
  # Ensure directory exists
  if (!dir.exists(DIR)) dir.create(DIR, recursive = TRUE)
  
  ggsave(
    filename = file.path(DIR, file_name),
    plot = plt,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300
  )
}

truncate_variants <- function(variants) {
  sub("^(?:[^:]*:){2}([^:]*):.*$", "\\1", variants)
}

