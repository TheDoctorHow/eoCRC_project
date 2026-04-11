# mechanism/scripts/fetch_string_network.R
#
# Fetch STRING PPI subnetworks for human and bacterial seed proteins.
# Reads: mechanism/data/seed_proteins.tsv
# Writes per-organism edge + node TSVs to mechanism/results/
#
# Run from repo root:
#   Rscript mechanism/scripts/fetch_string_network.R

# ---------------------------------------------------------------------------
# 0. Install dependencies if absent
# ---------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("STRINGdb", quietly = TRUE))
  BiocManager::install("STRINGdb")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

suppressPackageStartupMessages({
  library(STRINGdb)
  library(dplyr)
})

# ---------------------------------------------------------------------------
# 1. Paths — resolve repo root regardless of invocation style
# ---------------------------------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
flag_idx    <- grep("^--file=", args)
repo_root   <- if (length(flag_idx) == 1L) {
  script_path <- sub("^--file=", "", args[flag_idx])
  normalizePath(file.path(dirname(script_path), "../.."))
} else {
  getwd()   # fallback: interactive session, assume working dir = repo root
}

seed_file   <- file.path(repo_root, "mechanism", "data", "seed_proteins.tsv")
results_dir <- file.path(repo_root, "mechanism", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 2. Configuration
# ---------------------------------------------------------------------------
STRING_VERSION       <- "11.5"
STRING_SCORE_THRESH  <- 700L
MAX_SHELL_NODES      <- 10L   # max new (non-seed) nodes added per organism

# NCBI taxonomy IDs for each organism in the seed list
TAX_IDS <- c(
  "Homo sapiens"           = 9606L,
  "Parvimonas micra"       = 33033L,
  "Parvimonas stomatis"    = 341694L,
  "Gemella morbillorum"    = 29390L,
  "Dialister pneumosintes" = 39950L
)

# Output filename slug (human gets the literal label "human"; bacteria get genus_species)
FILE_SLUG <- c(
  "Homo sapiens"           = "human",
  "Parvimonas micra"       = "parvimonas_micra",
  "Parvimonas stomatis"    = "parvimonas_stomatis",
  "Gemella morbillorum"    = "gemella_morbillorum",
  "Dialister pneumosintes" = "dialister_pneumosintes"
)

# ---------------------------------------------------------------------------
# 3. Read seed list (comment lines starting with # are stripped automatically)
# ---------------------------------------------------------------------------
seeds_raw <- read.table(
  seed_file,
  sep              = "\t",
  header           = TRUE,
  comment.char     = "#",
  stringsAsFactors = FALSE,
  quote            = ""
)

seeds <- seeds_raw[seeds_raw$organism %in% names(TAX_IDS), ]
cat(sprintf("Loaded %d seed entries across %d organisms.\n\n",
            nrow(seeds), length(unique(seeds$organism))))

# ---------------------------------------------------------------------------
# 4. Core fetcher: one organism's subnetwork + one shell of ≤ MAX_SHELL_NODES
# ---------------------------------------------------------------------------
fetch_organism_network <- function(org_name, tax_id, seed_df) {

  cat(sprintf("--- %s (taxid %d) ---\n", org_name, tax_id))

  proteins   <- seed_df$protein
  uniprots   <- seed_df$uniprot_id  # may contain NA

  # Warn about (and note) NA-UniProt seeds up front; we still attempt name-based
  # mapping for them because STRINGdb maps by gene symbol, not UniProt
  na_mask <- is.na(uniprots) | uniprots == "NA"
  if (any(na_mask)) {
    warning(sprintf(
      "[%s] %d seed protein(s) have NA UniProt ID — will attempt gene-symbol mapping anyway, but may fail: %s",
      org_name, sum(na_mask),
      paste(proteins[na_mask], collapse = ", ")
    ))
  }

  # -- Initialise STRINGdb (downloads species mapping on first call; cached) --
  cache_dir <- file.path(results_dir, sprintf("string_cache_%d", tax_id))
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  sdb <- tryCatch(
    STRINGdb$new(
      version         = STRING_VERSION,
      species         = tax_id,
      score_threshold = STRING_SCORE_THRESH,
      input_directory = cache_dir
    ),
    error = function(e) {
      warning(sprintf("[%s] STRINGdb init failed: %s", org_name, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(sdb)) {
    return(list(edges = NULL, nodes = NULL,
                failed = proteins, n_mapped = 0L))
  }

  # -- Map seed proteins by gene/protein name --------------------------------
  query_df <- data.frame(protein = proteins, stringsAsFactors = FALSE)

  mapped <- tryCatch(
    sdb$map(query_df, "protein", removeUnmappedRows = TRUE),
    error = function(e) {
      warning(sprintf("[%s] Mapping failed: %s", org_name, conditionMessage(e)))
      NULL
    }
  )

  # Normalise: older STRINGdb versions use column "STRING_id"; newer may differ
  if (!is.null(mapped) && !"STRING_id" %in% colnames(mapped) &&
      "string_id" %in% tolower(colnames(mapped))) {
    colnames(mapped)[tolower(colnames(mapped)) == "string_id"] <- "STRING_id"
  }

  if (is.null(mapped) || nrow(mapped) == 0L) {
    warning(sprintf("[%s] No seed proteins mapped to STRING IDs. Skipping.", org_name))
    return(list(edges = NULL, nodes = NULL,
                failed = proteins, n_mapped = 0L))
  }

  seed_ids     <- unique(mapped$STRING_id)
  mapped_names <- unique(mapped$protein)
  failed_names <- setdiff(proteins, mapped_names)

  cat(sprintf("  Mapped %d / %d seeds to STRING IDs.\n",
              length(seed_ids), length(proteins)))
  if (length(failed_names) > 0)
    cat(sprintf("  Unmapped: %s\n", paste(failed_names, collapse = ", ")))

  # -- Seed-level interaction edges ------------------------------------------
  seed_edges <- tryCatch(
    sdb$get_interactions(seed_ids),
    error = function(e) {
      warning(sprintf("[%s] get_interactions failed: %s", org_name, conditionMessage(e)))
      data.frame()
    }
  )

  # -- One-shell expansion ---------------------------------------------------
  shell_ids <- character(0)

  raw_neighbors <- tryCatch(
    sdb$get_neighbors(seed_ids),
    error = function(e) {
      warning(sprintf("[%s] get_neighbors failed: %s", org_name, conditionMessage(e)))
      character(0)
    }
  )

  new_neighbors <- setdiff(raw_neighbors, seed_ids)

  if (length(new_neighbors) > 0L) {
    # Score each new neighbor: max combined_score with any seed
    expanded_edges <- tryCatch(
      sdb$get_interactions(c(seed_ids, new_neighbors)),
      error = function(e) data.frame()
    )

    if (nrow(expanded_edges) > 0L) {
      neighbor_scores <- expanded_edges %>%
        filter(
          (from %in% new_neighbors & to   %in% seed_ids) |
          (to   %in% new_neighbors & from %in% seed_ids)
        ) %>%
        mutate(neighbor_id = if_else(from %in% new_neighbors, from, to)) %>%
        group_by(neighbor_id) %>%
        summarise(max_score = max(combined_score, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_score))

      shell_ids <- head(neighbor_scores$neighbor_id, MAX_SHELL_NODES)
      cat(sprintf("  Shell expansion: %d new nodes added (cap = %d).\n",
                  length(shell_ids), MAX_SHELL_NODES))
    }
  }

  # -- Full edge + node tables for seeds ∪ shell -----------------------------
  all_ids <- unique(c(seed_ids, shell_ids))

  full_edges <- if (length(shell_ids) > 0L) {
    tryCatch(
      sdb$get_interactions(all_ids),
      error = function(e) { warning(conditionMessage(e)); data.frame() }
    )
  } else {
    seed_edges
  }

  # Node table: join STRING protein info
  string_proteins <- tryCatch(
    sdb$get_proteins(),
    error = function(e) data.frame()
  )

  node_tbl <- data.frame(STRING_id = all_ids, stringsAsFactors = FALSE)

  if (nrow(string_proteins) > 0L) {
    # Normalise column name
    if (!"STRING_id" %in% colnames(string_proteins) &&
        "string_id" %in% tolower(colnames(string_proteins))) {
      colnames(string_proteins)[tolower(colnames(string_proteins)) == "string_id"] <- "STRING_id"
    }
    if ("STRING_id" %in% colnames(string_proteins)) {
      node_tbl <- left_join(node_tbl, string_proteins, by = "STRING_id")
    }
  }

  # Annotate seed vs shell, and carry over original protein symbol for seeds
  node_tbl$node_type <- ifelse(node_tbl$STRING_id %in% seed_ids, "seed", "shell")
  node_tbl <- left_join(
    node_tbl,
    mapped[, c("protein", "STRING_id"), drop = FALSE],
    by = "STRING_id"
  )

  cat(sprintf("  Final subnetwork: %d nodes, %d edges.\n",
              nrow(node_tbl),
              if (is.null(full_edges) || nrow(full_edges) == 0L) 0L else nrow(full_edges)))

  list(
    edges    = full_edges,
    nodes    = node_tbl,
    failed   = failed_names,
    n_mapped = length(seed_ids)
  )
}

# ---------------------------------------------------------------------------
# 5. Main loop over organisms
# ---------------------------------------------------------------------------
summary_list <- list()

for (org in names(TAX_IDS)) {

  org_seeds <- seeds[seeds$organism == org, ]
  if (nrow(org_seeds) == 0L) {
    cat(sprintf("No seed entries for '%s' — skipping.\n\n", org))
    next
  }

  result <- fetch_organism_network(
    org_name = org,
    tax_id   = TAX_IDS[[org]],
    seed_df  = org_seeds
  )

  slug    <- FILE_SLUG[[org]]
  n_edges <- 0L
  n_nodes <- 0L

  if (!is.null(result)) {
    if (!is.null(result$edges) && nrow(result$edges) > 0L) {
      edge_path <- file.path(results_dir,
                             sprintf("string_edges_%s.tsv", slug))
      write.table(result$edges, edge_path,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      n_edges <- nrow(result$edges)
      cat(sprintf("  -> Saved %s\n", edge_path))
    } else {
      cat(sprintf("  No edges to save for %s.\n", org))
    }

    if (!is.null(result$nodes) && nrow(result$nodes) > 0L) {
      node_path <- file.path(results_dir,
                             sprintf("string_nodes_%s.tsv", slug))
      write.table(result$nodes, node_path,
                  sep = "\t", row.names = FALSE, quote = FALSE)
      n_nodes <- nrow(result$nodes)
      cat(sprintf("  -> Saved %s\n", node_path))
    }

    summary_list[[org]] <- data.frame(
      organism       = org,
      tax_id         = TAX_IDS[[org]],
      seeds_input    = nrow(org_seeds),
      seeds_mapped   = result$n_mapped,
      n_nodes        = n_nodes,
      n_edges        = n_edges,
      failed_mapping = if (length(result$failed) > 0L)
                         paste(result$failed, collapse = "; ")
                       else "none",
      stringsAsFactors = FALSE
    )
  }
  cat("\n")
}

# ---------------------------------------------------------------------------
# 6. Print and save summary
# ---------------------------------------------------------------------------
cat("==================== SUMMARY ====================\n")
if (length(summary_list) > 0L) {
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  print(summary_df, right = FALSE)

  summary_path <- file.path(results_dir, "string_fetch_summary.tsv")
  write.table(summary_df, summary_path,
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("\nSummary written to %s\n", summary_path))
} else {
  cat("No results produced — check warnings above.\n")
}
cat("=================================================\n")
