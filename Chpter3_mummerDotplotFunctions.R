# Functions for reading and plotting MUMmer delta alignments.

# Read a delta file and return alignment coordinates.
readDelta <- function(deltafile) {
  lines <- scan(deltafile, "a", sep = "\n", quiet = TRUE)
  lines <- lines[-1]
  lines_l <- strsplit(lines, " ")
  lines_len <- lapply(lines_l, length) %>% as.numeric()
  lines_l <- lines_l[lines_len != 1]
  lines_len <- lines_len[lines_len != 1]
  head_pos <- which(lines_len == 4)
  head_id <- rep(head_pos, c(head_pos[-1], length(lines_l) + 1) - head_pos)
  mat <- matrix(as.numeric(unlist(lines_l[lines_len == 7])), 7)
  res <- as.data.frame(t(mat[1:5, ]))
  colnames(res) <- c("rs", "re", "qs", "qe", "error")
  res$qid <- unlist(lapply(lines_l[head_id[lines_len == 7]], "[", 2))
  res$rid <- unlist(lapply(lines_l[head_id[lines_len == 7]], "[", 1)) %>%
    gsub("^>", "", .)
  res$strand <- ifelse(res$qe - res$qs > 0, "+", "-")
  res
}

# Calculate alignment length, percentage identity and midpoint.
calculate_perc_id_mid_points <- function(data) {
  data$r_length <- data$re - data$rs
  data$perc_id <- ((data$r_length - data$error) / data$r_length) * 100
  data$perc_id_factor <- data$perc_id
  data[data$perc_id < 100, "perc_id_factor"] <- "<100"
  data$perc_id_factor <- as.factor(data$perc_id_factor)
  data$r_mid <- (data$rs + data$re) / 2
  data$q_mid <- (data$qs + data$qe) / 2
  data
}

# Read a delta file and filter by minimum alignment length.
pre_plot_analysis <- function(delta_path, min_size = 20000) {
  data <- readDelta(delta_path)
  data <- calculate_perc_id_mid_points(data)
  data[data$r_length >= min_size, ]
}

# Read a delta file and filter by length and reference interval.
pre_plot_analysis_ref_filt <- function(delta_path, min_size = 20000,
                                       ref_start, ref_end) {
  data <- pre_plot_analysis(delta_path, min_size = min_size)
  data[data$rs >= ref_start & data$re <= ref_end, ]
}

# Plot pairwise alignments coloured by percentage identity.
plot_by_perc_id_cap <- function(data, xmin = 0, xmax = max(data$re),
                                cap_lower, cap_upper) {
  data <- data[order(data$perc_id), ]

  ggplot(
    data[data$rs > xmin & data$re < xmax, ],
    aes(x = rs, xend = re, y = qs, yend = qe, colour = perc_id)
  ) +
    geom_segment() +
    geom_point(alpha = 0.5, size = 0.4) +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 180, size = 5),
      strip.background = element_blank()
    ) +
    xlab(unique(data$rid)) +
    ylab(unique(data$qid)) +
    scale_colour_viridis(limits = c(cap_lower, cap_upper)) +
    scale_x_continuous(breaks = seq(xmin, xmax, 100000000)) +
    scale_y_continuous(breaks = seq(0, max(data$qe), 100000000))
}

# Plot alignment midpoints coloured by percentage identity.
plot_by_perc_id_cap_midpoint <- function(data, xmin = 0,
                                         xmax = max(data$re),
                                         cap_lower, cap_upper) {
  data <- data[order(data$perc_id), ]

  ggplot(
    data[data$rs > xmin & data$re < xmax, ],
    aes(x = r_mid, y = q_mid, colour = perc_id)
  ) +
    geom_point(size = 0.3) +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 180, size = 5),
      strip.background = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_colour_viridis(limits = c(cap_lower, cap_upper)) +
    scale_x_continuous(breaks = seq(xmin, xmax, 100000000)) +
    scale_y_continuous(breaks = seq(0, max(data$qe), 100000000)) +
    coord_cartesian(xlim = c(xmin, xmax))
}

# Plot percentage identity against the reference position.
plot_perc_id_v_ref <- function(data, xmin = 0, xmax = max(data$re),
                               ymin = 97, ymax = 100) {
  ggplot(
    data[data$r_mid > xmin & data$r_mid < xmax, ],
    aes(x = r_mid, y = perc_id, colour = r_length)
  ) +
    theme_bw() +
    xlab(data$rid) +
    ylab(paste0("percentage ID v ", data$qid)) +
    geom_point(alpha = 0.5) +
    ylim(ymin, ymax) +
    scale_colour_viridis()
}

# Plot percentage identity without colouring by alignment length.
plot_perc_id_v_ref_no_length <- function(data, xmin = 0,
                                         xmax = max(data$re),
                                         ymin = 97, ymax = 100,
                                         breaks = 0.5) {
  ggplot(
    data[data$r_mid > xmin & data$r_mid < xmax, ],
    aes(x = r_mid, y = perc_id)
  ) +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 180, size = 5),
      strip.background = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    geom_point(alpha = 0.5, size = 0.5, colour = "gray25") +
    scale_x_continuous(breaks = seq(xmin, xmax, 100000000)) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = seq(ymin, ymax, breaks)
    ) +
    coord_cartesian(xlim = c(xmin, xmax))
}

# Assign alignments to reference-position bins.
bin_data <- function(data, bin_size, max_chrom_size) {
  bins <- seq(0, max_chrom_size, by = bin_size)
  data$bin <- NA

  for (i in bins) {
    data$bin <- ifelse(
      data$r_mid > (i - bin_size) & data$r_mid < (i - 1),
      i,
      data$bin
    )
  }

  data
}

# Plot median percentage identity for each reference-position bin.
plot_medians_line_colour_cut_off <- function(data, bin_size = 10000000,
                                             cut_off = 99.99,
                                             ymin = 97, ymax = 100,
                                             max_chrom_size) {
  comparison_filt_bin <- bin_data(
    data,
    bin_size = bin_size,
    max_chrom_size = max_chrom_size
  )
  comparison_medians <- aggregate(
    perc_id ~ bin,
    data = comparison_filt_bin,
    FUN = median
  )
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- ifelse(
    comparison_medians$perc_id_median >= cut_off,
    paste0(">=", cut_off),
    paste0("<", cut_off)
  )
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)

  ggplot(
    comparison_medians,
    aes(x = bin, y = perc_id_median, colour = cut_off)
  ) +
    geom_line(colour = "grey") +
    geom_point(size = 1) +
    ylim(ymin, ymax) +
    scale_colour_manual(values = c("#73D055FF", "#440154FF")) +
    labs(colour = "% id") +
    theme_bw() +
    xlab(data$rid) +
    ggtitle(
      paste0(
        data$rid, " v ", data$qid,
        " BinSize: ", bin_size,
        " CutOff: ", cut_off
      )
    )
}

# Plot binned percentage identities as boxplots.
plot_boxplots_median_colour_cut_off <- function(data, bin_size = 10000000,
                                                cut_off = 99.99,
                                                ymin = 97, ymax = 100,
                                                max_chrom_size) {
  comparison_filt_bin <- bin_data(
    data,
    bin_size = bin_size,
    max_chrom_size = max_chrom_size
  )
  comparison_medians <- aggregate(
    perc_id ~ bin,
    data = comparison_filt_bin,
    FUN = median
  )
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- ifelse(
    comparison_medians$perc_id_median >= cut_off,
    paste0(">=", cut_off),
    paste0("<", cut_off)
  )
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)

  comparison_to_plot <- merge(comparison_filt_bin, comparison_medians)

  ggplot(
    comparison_to_plot,
    aes(x = bin, y = perc_id, group = bin, fill = cut_off)
  ) +
    geom_boxplot(outlier.shape = NA, lwd = 0.3) +
    scale_fill_manual(values = c("darkorchid4", "gold")) +
    theme_bw() +
    theme(
      strip.text.y = element_text(angle = 180, size = 5),
      strip.background = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    scale_x_continuous(breaks = seq(0, max(data$re), 100000000)) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = seq(ymin, ymax, 0.5)
    )
}
