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

pyseer_to_manhattan <- function(infile, chrom, outfile,
                                pval_col = "lrt-pvalue") {
  df <- read.delim(infile, header = TRUE, stringsAsFactors = FALSE)
  
  # extract position from variant name: CHROM_POS_REF_ALT -> POS
  pos <- as.numeric(sapply(strsplit(df$variant, "_"), `[`, 2))
  
  pval <- df[[pval_col]]
  neglog10p <- -log10(pval)
  
  out <- data.frame(
    `#CHR`       = chrom,
    SNP          = ".",
    BP           = pos,
    `minLOG10(P)`= neglog10p,
    `log10(p)`   = neglog10p,
    `r^2`        = 0,
    check.names = FALSE
  )
  
  # drop rows where p-value was NA/invalid (e.g. NaN, filtered variants)
  out <- out[is.finite(out$`minLOG10(P)`), ]
  
  write.table(out, outfile, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  invisible(out)
}

process_susie_results <- function(method_names, parameter_grid, iter,
                                  alpha = 0.75, data_dir) {
  
  n <- length(method_names) * nrow(parameter_grid) * iter
  results <- data.frame(
    method = character(n),
    iter = integer(n),
    grid_id = integer(n),
    fp = numeric(n),
    fdp = numeric(n),
    fwe = numeric(n),
    fdp_min = numeric(n),
    power_raw = numeric(n),
    power_bh = numeric(n),
    power_bonf = numeric(n),
    power_min = numeric(n),
    ssq_diff = numeric(n),
    tausq_diff = numeric(n)
  )
  
  count <- 1
  for (index in 1:nrow(parameter_grid)) {
    grid_id <- parameter_grid$grid_id[index]
    s <- parameter_grid$s[index]
    true_tausq <- parameter_grid$tausq[index]
    print(index)
    
    for (method in method_names) {
      for (j in 1:iter) {
        cur_result <- tryCatch(
          readRDS(paste0(data_dir, "/", method, "/", method, "_index_", index, "_iter_", j, "_results.rds")),
          error = function(e) NULL
        )
        
        if (is.null(cur_result) || (length(cur_result) == 1 && is.na(cur_result))) next
        
        nonnulls <- cur_result$name
        
        if (!is.list(cur_result)) {
          fp <- NA; fdp <- NA; fwe <- NA; fdp_min <- NA
          power_raw <- NA; power_bh <- NA; power_bonf <- NA; power_min <- NA
          ssq_diff <- NA; tausq_diff <- NA
        } else {
          if (is.null(cur_result$pip)) {
            p.vals <- 1 - cur_result$marginalPIP
            variants <- cur_result$variants
            est_ssq <- cur_result$sigmasq
            est_tausq <- cur_result$tausq
          } else {
            p.vals <- 1 - cur_result$pip
            variants <- names(cur_result$pip)
            est_ssq <- cur_result$sigma2
            est_tausq <- 0
          }
          
          nonnulls.raw <- variants[which(p.vals < alpha)]
          nonnulls.bh <- variants[which(p.adjust(p.vals, method = "BH") < alpha)]
          nonnulls.bonf <- variants[which(p.adjust(p.vals, method = "bonferroni") < alpha)]
          nonnulls.min <- variants[order(p.vals)[1:s]]
          
          fp <- length(setdiff(nonnulls.raw, nonnulls)) / max(1, length(nonnulls.raw))
          fdp <- length(setdiff(nonnulls.bh, nonnulls)) / max(1, length(nonnulls.bh))
          fwe <- length(setdiff(nonnulls.bonf, nonnulls)) >= 1
          fdp_min <- length(setdiff(nonnulls.min, nonnulls)) / max(1, length(nonnulls.min))
          power_raw <- sum(nonnulls %in% nonnulls.raw) / s
          power_bh <- sum(nonnulls %in% nonnulls.bh) / s
          power_bonf <- sum(nonnulls %in% nonnulls.bonf) / s
          power_min <- sum(nonnulls %in% nonnulls.min) / s
          ssq_diff <- est_ssq - 1
          tausq_diff <- est_tausq - true_tausq
        }
        
        results[count,] <- c(method, j, grid_id, fp, fdp, fwe, fdp_min,
                             power_raw, power_bh, power_bonf, power_min,
                             ssq_diff, tausq_diff)
        count <- count + 1
      }
    }
  }
  
  results$fp <- as.numeric(results$fp)
  results$fdp <- as.numeric(results$fdp)
  results$fwe <- as.logical(results$fwe)
  results$fdp_min <- as.numeric(results$fdp_min)
  results$power_raw <- as.numeric(results$power_raw)
  results$power_bh <- as.numeric(results$power_bh)
  results$power_bonf <- as.numeric(results$power_bonf)
  results$power_min <- as.numeric(results$power_min)
  results$ssq_diff <- as.numeric(results$ssq_diff)
  results$tausq_diff <- as.numeric(results$tausq_diff)
  results$grid_id <- as.integer(results$grid_id)
  results$iter <- as.integer(results$iter)
  
  results
}

process_pyseer_results <- function(method_names, parameter_grid, iter,
                                   alpha = 0.75, data_dir) {
  
  n <- length(method_names) * nrow(parameter_grid) * iter
  results <- data.frame(
    method = character(n),
    iter = integer(n),
    grid_id = integer(n),
    fp = numeric(n),
    fdp = numeric(n),
    fwe = numeric(n),
    fdp_min = numeric(n),
    power_raw = numeric(n),
    power_bh = numeric(n),
    power_bonf = numeric(n),
    power_min = numeric(n)
  )
  
  count <- 1
  for (index in 1:nrow(parameter_grid)) {
    grid_id <- parameter_grid$grid_id[index]
    s <- parameter_grid$s[index]
    print(index)
    
    for (method in method_names) {
      for (j in 1:iter) {
        cur_result <- tryCatch(
          readRDS(paste0(data_dir, "/", method, "/", method, "_index_", index, "_iter_", j, "_results.rds")),
          error = function(e) NULL
        )
        
        if (is.null(cur_result) || (length(cur_result) == 1 && is.na(cur_result))) next
        
        nonnulls <- cur_result$name
        
        has_pvals <- (is.data.frame(cur_result) || is.list(cur_result)) &&
          !is.null(cur_result$lrt.pvalue) &&
          !is.null(cur_result$variant)
        
        if (!has_pvals) {
          fp <- NA; fdp <- NA; fwe <- NA; fdp_min <- NA
          power_raw <- NA; power_bh <- NA; power_bonf <- NA; power_min <- NA
        } else {
          p.vals <- cur_result$lrt.pvalue
          variants <- cur_result$variant
          
          keep <- !is.na(p.vals)
          p.vals <- p.vals[keep]
          variants <- variants[keep]
          
          if (length(p.vals) == 0) {
            fp <- NA; fdp <- NA; fwe <- NA; fdp_min <- NA
            power_raw <- NA; power_bh <- NA; power_bonf <- NA; power_min <- NA
          } else {
            nonnulls.raw <- variants[which(p.vals < alpha)]
            nonnulls.bh <- variants[which(p.adjust(p.vals, method = "BH") < alpha)]
            nonnulls.bonf <- variants[which(p.adjust(p.vals, method = "bonferroni") < alpha)]
            nonnulls.min <- variants[order(p.vals)[1:s]]
            
            fp <- length(setdiff(nonnulls.raw, nonnulls)) / max(1, length(nonnulls.raw))
            fdp <- length(setdiff(nonnulls.bh, nonnulls)) / max(1, length(nonnulls.bh))
            fwe <- length(setdiff(nonnulls.bonf, nonnulls)) >= 1
            fdp_min <- length(setdiff(nonnulls.min, nonnulls)) / max(1, length(nonnulls.min))
            power_raw <- sum(nonnulls %in% nonnulls.raw) / s
            power_bh <- sum(nonnulls %in% nonnulls.bh) / s
            power_bonf <- sum(nonnulls %in% nonnulls.bonf) / s
            power_min <- sum(nonnulls %in% nonnulls.min) / s
          }
        }
        
        results[count,] <- c(method, j, grid_id, fp, fdp, fwe, fdp_min,
                             power_raw, power_bh, power_bonf, power_min)
        count <- count + 1
      }
    }
  }
  
  results$fp <- as.numeric(results$fp)
  results$fdp <- as.numeric(results$fdp)
  results$fwe <- as.logical(results$fwe)
  results$fdp_min <- as.numeric(results$fdp_min)
  results$power_raw <- as.numeric(results$power_raw)
  results$power_bh <- as.numeric(results$power_bh)
  results$power_bonf <- as.numeric(results$power_bonf)
  results$power_min <- as.numeric(results$power_min)
  results$grid_id <- as.integer(results$grid_id)
  results$iter <- as.integer(results$iter)
  
  results
}


process_pyseer_susie_results <- function(method_names, parameter_grid, iter,
                                         aggregation = c("mean", "median", "max"),
                                         alpha = 0.75, data_dir) {
  
  n <- length(method_names) * nrow(parameter_grid) * iter * length(aggregation)
  results <- data.frame(
    method = character(n),
    aggregation = character(n),
    iter = integer(n),
    grid_id = integer(n),
    fp = numeric(n),
    fdp = numeric(n),
    fwe = numeric(n),
    fdp_min = numeric(n),
    power_raw = numeric(n),
    power_bh = numeric(n),
    power_bonf = numeric(n),
    power_min = numeric(n)
  )
  
  count <- 1
  for (index in 1:nrow(parameter_grid)) {
    grid_id <- parameter_grid$grid_id[index]
    s <- parameter_grid$s[index]
    print(index)
    
    for (method in method_names) {
      for (j in 1:iter) {
        cur_result <- tryCatch(
          readRDS(paste0(data_dir, "/", method, "/", method, "_index_", index, "_iter_", j, "_results.rds")),
          error = function(e) NULL
        )
        
        if (is.null(cur_result) || (length(cur_result) == 1 && is.na(cur_result))) next
        
        nonnulls <- cur_result$name
        
        has_pip <- is.list(cur_result) &&
          !is.null(cur_result$pip) &&
          !is.null(cur_result$variant_name)
        
        if (has_pip) {
          raw_df <- data.frame(
            variant_name = cur_result$variant_name,
            pip = cur_result$pip
          )
          raw_df <- raw_df[!is.na(raw_df$pip), ]
        }
        
        for (agg in aggregation) {
          
          if (!has_pip || nrow(raw_df) == 0) {
            fp <- NA; fdp <- NA; fwe <- NA; fdp_min <- NA
            power_raw <- NA; power_bh <- NA; power_bonf <- NA; power_min <- NA
          } else {
            agg_fun <- match.fun(agg)
            
            agg_df <- raw_df %>%
              group_by(variant_name) %>%
              summarise(pip_agg = agg_fun(pip, na.rm = TRUE), .groups = "drop")
            
            p.vals <- 1 - agg_df$pip_agg
            variants <- agg_df$variant_name
            
            nonnulls.raw <- variants[which(p.vals < alpha)]
            nonnulls.bh <- variants[which(p.adjust(p.vals, method = "BH") < alpha)]
            nonnulls.bonf <- variants[which(p.adjust(p.vals, method = "bonferroni") < alpha)]
            nonnulls.min <- variants[order(p.vals)[1:s]]
            
            fp <- length(setdiff(nonnulls.raw, nonnulls)) / max(1, length(nonnulls.raw))
            fdp <- length(setdiff(nonnulls.bh, nonnulls)) / max(1, length(nonnulls.bh))
            fwe <- length(setdiff(nonnulls.bonf, nonnulls)) >= 1
            fdp_min <- length(setdiff(nonnulls.min, nonnulls)) / max(1, length(nonnulls.min))
            power_raw <- sum(nonnulls %in% nonnulls.raw) / s
            power_bh <- sum(nonnulls %in% nonnulls.bh) / s
            power_bonf <- sum(nonnulls %in% nonnulls.bonf) / s
            power_min <- sum(nonnulls %in% nonnulls.min) / s
          }
          
          results[count,] <- c(method, agg, j, grid_id, fp, fdp, fwe, fdp_min,
                               power_raw, power_bh, power_bonf, power_min)
          count <- count + 1
        }
      }
    }
  }
  
  results$fp <- as.numeric(results$fp)
  results$fdp <- as.numeric(results$fdp)
  results$fwe <- as.logical(results$fwe)
  results$fdp_min <- as.numeric(results$fdp_min)
  results$power_raw <- as.numeric(results$power_raw)
  results$power_bh <- as.numeric(results$power_bh)
  results$power_bonf <- as.numeric(results$power_bonf)
  results$power_min <- as.numeric(results$power_min)
  results$grid_id <- as.integer(results$grid_id)
  results$iter <- as.integer(results$iter)
  
  results
}

process_pyseer_susie_zero_results <- function(method_names, parameter_grid, iter,
                                              aggregation = c("mean", "median", "max"),
                                              alpha = 0.75, data_dir) {
  
  n <- length(method_names) * nrow(parameter_grid) * iter * length(aggregation)
  results <- data.frame(
    method = character(n),
    aggregation = character(n),
    iter = integer(n),
    grid_id = integer(n),
    fp = numeric(n),
    fdp = numeric(n),
    fwe = numeric(n),
    fdp_min = numeric(n),
    power_raw = numeric(n),
    power_bh = numeric(n),
    power_bonf = numeric(n),
    power_min = numeric(n)
  )
  
  count <- 1
  for (index in 1:nrow(parameter_grid)) {
    grid_id <- parameter_grid$grid_id[index]
    s <- parameter_grid$s[index]
    print(index)
    
    for (method in method_names) {
      for (j in 1:iter) {
        cur_result <- tryCatch(
          readRDS(paste0(data_dir, "/", method, "/", method, "_index_", index, "_iter_", j, "_results.rds")),
          error = function(e) NULL
        )
        
        if (is.null(cur_result) || (length(cur_result) == 1 && is.na(cur_result))) next
        
        nonnulls <- cur_result$name
        
        has_pip <- is.list(cur_result) &&
          !is.null(cur_result$pip) &&
          !is.null(cur_result$variant_name)
        
        if (has_pip) {
          raw_df <- data.frame(
            variant_name = cur_result$variant_name,
            pip = cur_result$pip
          )
          raw_df <- raw_df[!is.na(raw_df$pip), ]
        }
        
        for (agg in aggregation) {
          
          if (!has_pip || nrow(raw_df) == 0) {
            fp <- 0; fdp <- 0; fwe <- 0; fdp_min <- 0
            power_raw <- 0; power_bh <- 0; power_bonf <- 0; power_min <- 0
          } else {
            agg_fun <- match.fun(agg)
            
            agg_df <- raw_df %>%
              group_by(variant_name) %>%
              summarise(pip_agg = agg_fun(pip, na.rm = TRUE), .groups = "drop")
            
            p.vals <- 1 - agg_df$pip_agg
            variants <- agg_df$variant_name
            
            nonnulls.raw <- variants[which(p.vals < alpha)]
            nonnulls.bh <- variants[which(p.adjust(p.vals, method = "BH") < alpha)]
            nonnulls.bonf <- variants[which(p.adjust(p.vals, method = "bonferroni") < alpha)]
            nonnulls.min <- variants[order(p.vals)[1:s]]
            
            fp <- length(setdiff(nonnulls.raw, nonnulls)) / max(1, length(nonnulls.raw))
            fdp <- length(setdiff(nonnulls.bh, nonnulls)) / max(1, length(nonnulls.bh))
            fwe <- length(setdiff(nonnulls.bonf, nonnulls)) >= 1
            fdp_min <- length(setdiff(nonnulls.min, nonnulls)) / max(1, length(nonnulls.min))
            power_raw <- sum(nonnulls %in% nonnulls.raw) / s
            power_bh <- sum(nonnulls %in% nonnulls.bh) / s
            power_bonf <- sum(nonnulls %in% nonnulls.bonf) / s
            power_min <- sum(nonnulls %in% nonnulls.min) / s
          }
          
          results[count,] <- c(method, agg, j, grid_id, fp, fdp, fwe, fdp_min,
                               power_raw, power_bh, power_bonf, power_min)
          count <- count + 1
        }
      }
    }
  }
  
  results$fp <- as.numeric(results$fp)
  results$fdp <- as.numeric(results$fdp)
  results$fwe <- as.logical(results$fwe)
  results$fdp_min <- as.numeric(results$fdp_min)
  results$power_raw <- as.numeric(results$power_raw)
  results$power_bh <- as.numeric(results$power_bh)
  results$power_bonf <- as.numeric(results$power_bonf)
  results$power_min <- as.numeric(results$power_min)
  results$grid_id <- as.integer(results$grid_id)
  results$iter <- as.integer(results$iter)
  
  results
}

process_susie_zero_results <- function(method_names, parameter_grid, iter,
                                       alpha = 0.75, data_dir) {
  
  n <- length(method_names) * nrow(parameter_grid) * iter
  results <- data.frame(
    method = character(n),
    iter = integer(n),
    grid_id = integer(n),
    fp = numeric(n),
    fdp = numeric(n),
    fwe = numeric(n),
    fdp_min = numeric(n),
    power_raw = numeric(n),
    power_bh = numeric(n),
    power_bonf = numeric(n),
    power_min = numeric(n),
    ssq_diff = numeric(n),
    tausq_diff = numeric(n)
  )
  
  count <- 1
  for (index in 1:nrow(parameter_grid)) {
    grid_id <- parameter_grid$grid_id[index]
    s <- parameter_grid$s[index]
    true_tausq <- parameter_grid$tausq[index]
    print(index)
    
    for (method in method_names) {
      for (j in 1:iter) {
        cur_result <- tryCatch(
          readRDS(paste0(data_dir, "/", method, "/", method, "_index_", index, "_iter_", j, "_results.rds")),
          error = function(e) NULL
        )
        
        if (is.null(cur_result) || (length(cur_result) == 1 && is.na(cur_result))) next
        
        nonnulls <- cur_result$name
        
        if (!is.list(cur_result)) {
          fp <- 0; fdp <- 0; fwe <- 0; fdp_min <- 0
          power_raw <- 0; power_bh <- 0; power_bonf <- 0; power_min <- 0
          ssq_diff <- 0; tausq_diff <- 0
        } else {
          if (is.null(cur_result$pip)) {
            p.vals <- 1 - cur_result$marginalPIP
            variants <- cur_result$variants
            est_ssq <- cur_result$sigmasq
            est_tausq <- cur_result$tausq
          } else {
            p.vals <- 1 - cur_result$pip
            variants <- names(cur_result$pip)
            est_ssq <- cur_result$sigma2
            est_tausq <- 0
          }
          
          nonnulls.raw <- variants[which(p.vals < alpha)]
          nonnulls.bh <- variants[which(p.adjust(p.vals, method = "BH") < alpha)]
          nonnulls.bonf <- variants[which(p.adjust(p.vals, method = "bonferroni") < alpha)]
          nonnulls.min <- variants[order(p.vals)[1:s]]
          
          fp <- length(setdiff(nonnulls.raw, nonnulls)) / max(1, length(nonnulls.raw))
          fdp <- length(setdiff(nonnulls.bh, nonnulls)) / max(1, length(nonnulls.bh))
          fwe <- length(setdiff(nonnulls.bonf, nonnulls)) >= 1
          fdp_min <- length(setdiff(nonnulls.min, nonnulls)) / max(1, length(nonnulls.min))
          power_raw <- sum(nonnulls %in% nonnulls.raw) / s
          power_bh <- sum(nonnulls %in% nonnulls.bh) / s
          power_bonf <- sum(nonnulls %in% nonnulls.bonf) / s
          power_min <- sum(nonnulls %in% nonnulls.min) / s
          ssq_diff <- est_ssq - 1
          tausq_diff <- est_tausq - true_tausq
        }
        
        results[count,] <- c(method, j, grid_id, fp, fdp, fwe, fdp_min,
                             power_raw, power_bh, power_bonf, power_min,
                             ssq_diff, tausq_diff)
        count <- count + 1
      }
    }
  }
  
  results$fp <- as.numeric(results$fp)
  results$fdp <- as.numeric(results$fdp)
  results$fwe <- as.logical(results$fwe)
  results$fdp_min <- as.numeric(results$fdp_min)
  results$power_raw <- as.numeric(results$power_raw)
  results$power_bh <- as.numeric(results$power_bh)
  results$power_bonf <- as.numeric(results$power_bonf)
  results$power_min <- as.numeric(results$power_min)
  results$ssq_diff <- as.numeric(results$ssq_diff)
  results$tausq_diff <- as.numeric(results$tausq_diff)
  results$grid_id <- as.integer(results$grid_id)
  results$iter <- as.integer(results$iter)
  
  results
}

process_pyseer_zero_results <- function(method_names, parameter_grid, iter,
                                        alpha = 0.75, data_dir) {
  
  n <- length(method_names) * nrow(parameter_grid) * iter
  results <- data.frame(
    method = character(n),
    iter = integer(n),
    grid_id = integer(n),
    fp = numeric(n),
    fdp = numeric(n),
    fwe = numeric(n),
    fdp_min = numeric(n),
    power_raw = numeric(n),
    power_bh = numeric(n),
    power_bonf = numeric(n),
    power_min = numeric(n)
  )
  
  count <- 1
  for (index in 1:nrow(parameter_grid)) {
    grid_id <- parameter_grid$grid_id[index]
    s <- parameter_grid$s[index]
    print(index)
    
    for (method in method_names) {
      for (j in 1:iter) {
        cur_result <- tryCatch(
          readRDS(paste0(data_dir, "/", method, "/", method, "_index_", index, "_iter_", j, "_results.rds")),
          error = function(e) NULL
        )
        
        if (is.null(cur_result) || (length(cur_result) == 1 && is.na(cur_result))) next
        
        nonnulls <- cur_result$name
        
        has_pvals <- (is.data.frame(cur_result) || is.list(cur_result)) &&
          !is.null(cur_result$lrt.pvalue) &&
          !is.null(cur_result$variant)
        
        if (!has_pvals) {
          fp <- NA; fdp <- NA; fwe <- NA; fdp_min <- NA
          power_raw <- NA; power_bh <- NA; power_bonf <- NA; power_min <- NA
        } else {
          p.vals <- cur_result$lrt.pvalue
          variants <- cur_result$variant
          
          keep <- !is.na(p.vals)
          p.vals <- p.vals[keep]
          variants <- variants[keep]
          
          if (length(p.vals) == 0) {
            fp <- 0; fdp <- 0; fwe <- 0; fdp_min <- 0
            power_raw <- 0; power_bh <- 0; power_bonf <- 0; power_min <- 0
          } else {
            nonnulls.raw <- variants[which(p.vals < alpha)]
            nonnulls.bh <- variants[which(p.adjust(p.vals, method = "BH") < alpha)]
            nonnulls.bonf <- variants[which(p.adjust(p.vals, method = "bonferroni") < alpha)]
            nonnulls.min <- variants[order(p.vals)[1:s]]
            
            fp <- length(setdiff(nonnulls.raw, nonnulls)) / max(1, length(nonnulls.raw))
            fdp <- length(setdiff(nonnulls.bh, nonnulls)) / max(1, length(nonnulls.bh))
            fwe <- length(setdiff(nonnulls.bonf, nonnulls)) >= 1
            fdp_min <- length(setdiff(nonnulls.min, nonnulls)) / max(1, length(nonnulls.min))
            power_raw <- sum(nonnulls %in% nonnulls.raw) / s
            power_bh <- sum(nonnulls %in% nonnulls.bh) / s
            power_bonf <- sum(nonnulls %in% nonnulls.bonf) / s
            power_min <- sum(nonnulls %in% nonnulls.min) / s
          }
        }
        
        results[count,] <- c(method, j, grid_id, fp, fdp, fwe, fdp_min,
                             power_raw, power_bh, power_bonf, power_min)
        count <- count + 1
      }
    }
  }
  
  results$fp <- as.numeric(results$fp)
  results$fdp <- as.numeric(results$fdp)
  results$fwe <- as.logical(results$fwe)
  results$fdp_min <- as.numeric(results$fdp_min)
  results$power_raw <- as.numeric(results$power_raw)
  results$power_bh <- as.numeric(results$power_bh)
  results$power_bonf <- as.numeric(results$power_bonf)
  results$power_min <- as.numeric(results$power_min)
  results$grid_id <- as.integer(results$grid_id)
  results$iter <- as.integer(results$iter)
  
  results
}