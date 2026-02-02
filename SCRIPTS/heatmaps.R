#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

# =========================================================
# USAGE (flexible)
# ---------------------------------------------------------
# Minimal (no map / no CDS/mRNA):
#   Rscript heeatmaps.R <input.tsv> <output_dir>
#
# With optional mapping (SRR -> cell line):
#   Rscript heeatmaps.R <input.tsv> <output_dir> --map <map.tsv>
#
# With optional CDS/mRNA FASTA for start/stop vlines:
#   Rscript heeatmaps.R <input.tsv> <output_dir> --cds <cds.fasta> --mrna <mrna.fasta>
#
# Optional:
#   --offset <int>
#   --exclude "SRRxxx_1,SRRyyy_2"  (suffix _<n> ignored for matching)
#
# NOTES:
# - If --map is missing, heatmaps use the raw SRR/ERR column names (cleaned).
# - If --cds/--mrna are missing, no start/stop lines are drawn.
# =========================================================

# --------- X-axis label settings ---------
TOP_N_LABELS <- 2L               # number of density peaks to label (0 to disable)
MIN_GAP_BETWEEN_LABELS <- 15L    # min distance between labeled peaks (in k-mer indices)
EVERY_N_TICKS <- 0L              # 0 = off; otherwise adds a tick every N k-mers

KMER_LEN <- 25L
STRIDE   <- 1L

# ---------- CLI parsing ----------
args <- commandArgs(trailingOnly = TRUE)

stop_usage <- function() {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript heeatmaps.R <input.tsv> <output_dir> [--map map.tsv] [--cds cds.fasta --mrna mrna.fasta] [--offset int] [--exclude list]\n",
      "Examples:\n",
      "  Rscript heeatmaps.R input.tsv plots/\n",
      "  Rscript heeatmaps.R input.tsv plots/ --map map.tsv\n",
      "  Rscript heeatmaps.R input.tsv plots/ --cds cds.fa --mrna mrna.fa --offset 1\n"
    ),
    call. = FALSE
  )
}

if (length(args) < 2) stop_usage()

infile <- args[1]
outdir <- args[2]

# defaults
mapfile <- NA_character_
cds_fa  <- NA_character_
mrna_fa <- NA_character_
OFFSET  <- 0L
EXCLUDE_ARG <- ""

# parse flags
i <- 3
while (i <= length(args)) {
  a <- args[i]
  if (a == "--map") {
    i <- i + 1
    if (i > length(args)) stop("Missing value after --map", call. = FALSE)
    mapfile <- args[i]
  } else if (a == "--cds") {
    i <- i + 1
    if (i > length(args)) stop("Missing value after --cds", call. = FALSE)
    cds_fa <- args[i]
  } else if (a == "--mrna") {
    i <- i + 1
    if (i > length(args)) stop("Missing value after --mrna", call. = FALSE)
    mrna_fa <- args[i]
  } else if (a == "--offset") {
    i <- i + 1
    if (i > length(args)) stop("Missing value after --offset", call. = FALSE)
    OFFSET <- suppressWarnings(as.integer(args[i]))
    if (is.na(OFFSET)) OFFSET <- 0L
  } else if (a == "--exclude") {
    i <- i + 1
    if (i > length(args)) stop("Missing value after --exclude", call. = FALSE)
    EXCLUDE_ARG <- args[i]
  } else {
    stop("Unknown argument: ", a, call. = FALSE)
  }
  i <- i + 1
}

if (!file.exists(infile)) stop("Input file not found: ", infile, call. = FALSE)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
strip_suffix <- function(x) sub("(_\\d+)$", "", x)   # SRRxxxxx_1 -> SRRxxxxx

read_fasta_map <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(list())
  txt <- readLines(path, warn = FALSE)
  if (!length(txt)) return(list())
  idx <- which(startsWith(txt, ">")); idx <- c(idx, length(txt) + 1)
  out <- list()
  for (j in seq_len(length(idx) - 1)) {
    header <- sub("^>", "", txt[idx[j]])
    name <- strsplit(header, "\\s+")[[1]][1]
    s <- paste(txt[(idx[j]+1):(idx[j+1]-1)], collapse = "")
    s <- toupper(gsub("\\s+", "", s))
    s <- chartr("U", "T", s); s <- gsub("[^ACGT]", "", s)
    out[[name]] <- s
  }
  out
}

locate_cds_on_mrna <- function(cds, mrna) {
  if (is.null(cds) || is.null(mrna) || nchar(cds) < 3 || nchar(mrna) < 3) return(NULL)
  m <- regexpr(cds, mrna, fixed = TRUE)
  if (m[1] == -1) return(NULL)
  start_mrna <- as.integer(m[1])
  L <- nchar(cds); last3 <- if (L >= 3) substr(cds, L-2, L) else ""
  if (last3 %in% c("TAA","TAG","TGA")) {
    list(start=start_mrna, stop=start_mrna+L-3, stop_codon=last3, stop_included=TRUE)
  } else {
    list(start=start_mrna, stop=start_mrna+L-1, stop_codon=NA_character_, stop_included=FALSE)
  }
}

pos_to_kidx <- function(pos_nt, stride = STRIDE) floor((pos_nt - 1L)/stride) + 1L

# ---------- optional mapping SRR/ERR -> cell line ----------
use_map <- (!is.na(mapfile) && nzchar(mapfile) && file.exists(mapfile))
map_df <- NULL
if (use_map) {
  map_df <- read.table(mapfile, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("sample_id", "cell_line"))
  map_df$sample_id <- strip_suffix(map_df$sample_id)
} else {
  message("[INFO] No --map provided (or file not found). Using SRR/ERR names directly.")
}

# ---------- table kmers ----------
df <- read.csv(infile, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("id","tag") %in% names(df))) stop("Columns 'id' and 'tag' are required.", call. = FALSE)

# detect sample columns
sum_col <- which(names(df) == "sum")
if (length(sum_col) == 0) sum_col <- ncol(df) + 1
pat_cols <- 3:(sum_col - 1)
if (length(pat_cols) < 1) stop("No sample columns detected between columns 3 and 'sum'.", call. = FALSE)

# parse exclude list
excluded_raw <- character(0)
if (nzchar(EXCLUDE_ARG)) {
  excluded_raw <- unlist(strsplit(EXCLUDE_ARG, "[,;]"))
  excluded_raw <- trimws(excluded_raw)
  excluded_raw <- excluded_raw[nzchar(excluded_raw)]
}
excluded_clean <- strip_suffix(excluded_raw)

all_pat_names <- names(df)[pat_cols]
all_pat_clean <- strip_suffix(all_pat_names)
drop_idx <- which(all_pat_names %in% excluded_raw | all_pat_clean %in% excluded_clean)

if (length(drop_idx)) {
  message("[INFO] Excluding sample columns: ", paste(all_pat_names[drop_idx], collapse = ", "))
  keep <- setdiff(seq_along(all_pat_names), drop_idx)
  if (!length(keep)) stop("After exclusion, no sample columns remain.", call. = FALSE)
  pat_cols <- 2 + keep  # re-map to df indices (since pat_cols starts at 3)
}

# cast numeric
df[pat_cols] <- lapply(df[pat_cols], function(x) suppressWarnings(as.numeric(x)))

df <- df %>%
  mutate(
    contig = sub("_kmer_.*$", "", id),
    k_idx  = suppressWarnings(as.integer(sub("^.*_kmer_", "", id))),
    tag    = toupper(tag)
  )

# ---------- start/stop vlines (optional) ----------
cds_map  <- read_fasta_map(cds_fa)
mrna_map <- read_fasta_map(mrna_fa)

build_vlines <- function(contigs, df_kmers, cds_map, mrna_map, K = KMER_LEN) {
  if (!length(cds_map) || !length(mrna_map)) {
    return(data.frame(contig=character(0), start_k=integer(0), stop_k=integer(0), stringsAsFactors = FALSE))
  }
  recs <- lapply(contigs, function(ctg) {
    if (!ctg %in% names(cds_map) || !ctg %in% names(mrna_map)) return(NULL)
    cds  <- cds_map[[ctg]]; mrna <- mrna_map[[ctg]]
    if (nchar(cds) < K) return(NULL)

    start25 <- substr(cds, 1, K)
    dctg <- df_kmers %>% dplyr::filter(contig == ctg)
    hit_start <- dctg %>% dplyr::filter(tag == start25) %>% dplyr::arrange(k_idx)
    start_k <- if (nrow(hit_start) >= 1) hit_start$k_idx[1] else NA_integer_

    loc <- locate_cds_on_mrna(cds, mrna)
    if (is.null(loc)) {
      return(data.frame(contig=ctg, start_k=start_k, stop_k=NA_integer_, stringsAsFactors=FALSE))
    }

    stop_k <- NA_integer_
    if (loc$stop_included) {
      if (loc$stop + K - 1 <= nchar(mrna)) {
        stop25 <- substr(mrna, loc$stop, loc$stop + K - 1)
        hit_stop <- dctg %>% dplyr::filter(tag == stop25) %>% dplyr::arrange(k_idx)
        if (nrow(hit_stop) >= 1) stop_k <- hit_stop$k_idx[1]
      }
      if (is.na(stop_k)) stop_k <- pos_to_kidx(loc$stop)
    } else {
      endpos <- loc$stop
      if (endpos + K - 1 <= nchar(mrna)) {
        stop25 <- substr(mrna, endpos, endpos + K - 1)
        hit_stop <- dctg %>% dplyr::filter(tag == stop25) %>% dplyr::arrange(k_idx)
        if (nrow(hit_stop) >= 1) stop_k <- hit_stop$k_idx[1]
      }
      if (is.na(stop_k)) stop_k <- pos_to_kidx(endpos)
    }

    data.frame(contig=ctg, start_k=start_k, stop_k=stop_k, stringsAsFactors=FALSE)
  })
  dplyr::bind_rows(recs)
}

vlines <- build_vlines(unique(df$contig), df, cds_map, mrna_map, K = KMER_LEN)
if (!nrow(vlines)) message("[INFO] No CDS/mRNA provided (or not found). Skipping start/stop vlines.")

# ---------- long + z-score + labeling ----------
long <- df %>%
  select(contig, k_idx, tag, all_of(names(df)[pat_cols])) %>%
  pivot_longer(cols = -c(contig, k_idx, tag), names_to = "patient", values_to = "count") %>%
  mutate(patient_clean = strip_suffix(patient)) %>%
  {
    if (use_map) {
      left_join(., map_df, by = c("patient_clean" = "sample_id")) %>%
        mutate(patient_lab = ifelse(is.na(cell_line), patient_clean, cell_line))
    } else {
      mutate(., patient_lab = patient_clean)
    }
  } %>%
  arrange(contig, k_idx) %>%
  group_by(contig, patient_lab) %>%
  mutate(
    mu = mean(count, na.rm = TRUE),
    sdv = sd(count, na.rm = TRUE),
    z   = ifelse(is.finite(sdv) & sdv > 0, (count - mu)/sdv, 0)
  ) %>%
  ungroup() %>%
  mutate(
    z = pmax(-3, pmin(3, z)),
    col_group = ((k_idx - 1L) %% 3L) + 1L,
    grp = factor(col_group, levels = c(1,2,3), labels = c("p0","p1","p2")),
    a   = ifelse(z > 0, pmin(1, pmax(0, z/3)), 0),
    x_pos = k_idx + OFFSET
  )

# keep your colors
grp_cols <- c("p0"="#E66100","p1"="#005AC8","p2"="#C20088")

find_top_peaks <- function(d, n = 2L, min_gap = 15L) {
  if (n <= 0 || nrow(d) < 3) return(integer())
  y <- d$density
  locmax <- which(y > dplyr::lag(y, default = -Inf) & y >= dplyr::lead(y, default = -Inf))
  if (!length(locmax)) return(integer())
  cand <- d[locmax, , drop = FALSE] %>% arrange(desc(density))
  picks <- integer()
  for (ii in seq_len(nrow(cand))) {
    xi <- cand$x_pos[ii]
    if (!length(picks) || min(abs(xi - picks)) >= min_gap) {
      picks <- c(picks, xi)
      if (length(picks) >= n) break
    }
  }
  picks
}

# ---------- plot per contig ----------
plot_one <- function(df_contig, outdir) {
  ctg <- unique(df_contig$contig)
  n_p <- dplyr::n_distinct(df_contig$patient_lab)

  width_in  <- max(6, min(20, 0.22 * dplyr::n_distinct(df_contig$k_idx) + 3))
  height_in <- max(4.8, min(12.5, 0.26 * n_p + 1.4))

  dens_df <- df_contig %>%
    mutate(pos = pmax(z, 0)) %>%
    group_by(patient_lab) %>%
    mutate(
      scale_fac = suppressWarnings(quantile(pos[pos > 0], 0.95, na.rm = TRUE)),
      scale_fac = ifelse(is.finite(scale_fac) & scale_fac > 0, scale_fac, 1),
      pos_norm  = pmin(1, pos / scale_fac)
    ) %>%
    ungroup() %>%
    group_by(x_pos) %>%
    summarise(density = mean(pos_norm, na.rm = TRUE), .groups = "drop")

  map_ids <- df_contig %>% distinct(k_idx, x_pos) %>%
    mutate(label = paste0(ctg, "_kmer_", k_idx, " (x=", x_pos, ")"))

  peak_x <- find_top_peaks(dens_df, n = TOP_N_LABELS, min_gap = MIN_GAP_BETWEEN_LABELS)

  vv <- vlines[vlines$contig == ctg, , drop = FALSE]
  start_x <- if (nrow(vv) && !is.na(vv$start_k)) vv$start_k else integer()
  stop_x  <- if (nrow(vv) && !is.na(vv$stop_k))  vv$stop_k  else integer()

  reg_x <- integer()
  if (EVERY_N_TICKS > 0L) {
    rng <- range(df_contig$x_pos, na.rm = TRUE)
    reg_x <- seq(from = rng[1], to = rng[2], by = EVERY_N_TICKS)
  }

  breaks_x <- sort(unique(c(start_x, stop_x, peak_x, reg_x)))
  lab_tbl <- map_ids %>% filter(x_pos %in% breaks_x)
  label_x <- vapply(breaks_x, function(xx) {
    if (xx %in% c(start_x, stop_x, peak_x)) {
      lab_tbl$label[match(xx, lab_tbl$x_pos)]
    } else {
      as.character(xx)
    }
  }, character(1))

  p_top <- ggplot(dens_df, aes(x = x_pos, y = density)) +
    geom_area(alpha = 0.25) +
    geom_line(linewidth = 0.4) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                       breaks = breaks_x, labels = label_x) +
    labs(y = "density (norm.)", x = NULL, title = ctg) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title  = element_text(size = 13, face = "bold", margin = margin(b = 4)),
      axis.title.y= element_text(size = 9),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_blank(),
      panel.grid  = element_blank(),
      plot.margin = margin(6, 10, 0, 8)
    )

  if (length(start_x)) p_top <- p_top + geom_vline(xintercept = start_x, colour = "black", linewidth = 0.3, alpha = 0.35)
  if (length(stop_x))  p_top <- p_top  + geom_vline(xintercept = stop_x,  colour = "black", linewidth = 0.3, linetype = "dashed", alpha = 0.35)

  p_heat <- ggplot(df_contig, aes(x = x_pos, y = patient_lab)) +
    geom_tile(aes(fill = grp, alpha = a), width = 0.95, height = 0.95) +
    scale_fill_manual(values = grp_cols, name = "Phase") +
    scale_alpha(range = c(0, 1), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                       breaks = breaks_x, labels = label_x) +
    labs(
      x = if (OFFSET == 0L) "k-mer index" else paste0("k-mer index (offset=", OFFSET, ")"),
      y = if (use_map) "Cell line" else "Sample"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      panel.grid   = element_blank(),
      plot.margin  = margin(0, 10, 8, 8),
      legend.position = "right"
    ) +
    coord_cartesian(clip = "off")

  if (length(start_x)) p_heat <- p_heat + geom_vline(xintercept = start_x, colour = "black", linewidth = 0.3, alpha = 0.35, lineend = "round")
  if (length(stop_x))  p_heat <- p_heat  + geom_vline(xintercept = stop_x,  colour = "black", linewidth = 0.3, linetype = "dashed", alpha = 0.35, lineend = "round")

  p_final <- p_top / p_heat + patchwork::plot_layout(heights = c(1, 7))

  outfile <- file.path(outdir, paste0(ctg, "_heatmap_plus_density_offset", OFFSET, ".pdf"))
  pdf(outfile, width = width_in, height = height_in + 1.6)
  print(p_final)
  dev.off()
  message("✅ ", outfile)
}

invisible(lapply(split(long, long$contig), plot_one, outdir = outdir))
