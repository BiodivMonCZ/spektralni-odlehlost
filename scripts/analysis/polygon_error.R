# ==============================================================================
# polygon_error.R
# ==============================================================================

library(dplyr)
library(robust) 
library(bigutilsr)
library(sf)
library(purrr)

# 1. Nastavení relativních cest
dir_interim <- "interim"
dir_outputs <- "outputs"

file_stats <- file.path(dir_interim, "VMB_stats.csv")
file_shp   <- file.path(dir_interim, "VMB_filt.shp")
out_shp    <- file.path(dir_outputs, "VMB_outliers_polygon.shp")

if (!file.exists(file_stats)) stop("Chybí vstupní soubor: ", file_stats)
if (!file.exists(file_shp)) stop("Chybí vstupní soubor: ", file_shp)

# 2. Načtení a předzpracování dat
cat("Načítám tabulku statistik...\n")
df <- read.csv(file_stats)

df_model <- df %>%
  select(SEGMENT_ID, BIOTOP_SEZ, matches("mean|median|sd", ignore.case = TRUE)) %>%
  mutate(SEGMENT_ID = as.character(SEGMENT_ID),
         outlier = FALSE)

# 3. Identifikace validních biotopů (min 31 záznamů)
biotop_counts <- df_model %>% 
  count(BIOTOP_SEZ, name = "count") %>% 
  filter(count >= 31)

valid_biotops <- biotop_counts$BIOTOP_SEZ
cat("Nalezeno", length(valid_biotops), "kategorií biotopů splňujících limit (>= 31 polygonů).\n")

# 4. Výpočet odlehlostí
cat("Spouštím detekci odlehlostí (PCA + D2 + LOF)...\n")

detect_outliers <- function(biotop_name) {
  d <- df_model %>% filter(BIOTOP_SEZ == biotop_name)
  
  # Matice prediktorů
  X <- d %>% select(-SEGMENT_ID, -outlier, -BIOTOP_SEZ)
  
  # Odstranění sloupců s nulovým rozptylem
  X_var <- apply(X, 2, var, na.rm = TRUE)
  X_clean <- X[, X_var > 1e-6, drop = FALSE]
  
  if (ncol(X_clean) < 2) return(NULL)
  
  # Ošetření řádků bez NA v aktivních prediktorech
  complete_idx <- complete.cases(X_clean)
  if (sum(complete_idx) < 31) return(NULL)
  
  d_clean <- d[complete_idx, , drop = FALSE]
  X_clean <- X_clean[complete_idx, , drop = FALSE]
  
  # PCA
  pca <- tryCatch(prcomp(X_clean, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  
  n_comp <- min(10, ncol(pca$x))
  m <- pca$x[, 1:n_comp, drop = FALSE]
  
  # Výpočet původních metrik: D2 (dist_ogk) a LOF
  mah_dist <- tryCatch(bigutilsr::dist_ogk(m), error = function(e) rep(NA_real_, nrow(m)))
  lof_dist <- tryCatch(bigutilsr::LOF(m),      error = function(e) rep(NA_real_, nrow(m)))
  
  cand_MAH <- c(
    tryCatch(bigutilsr::hist_out(mah_dist)$lim[2], error = function(e) NA_real_),
    tryCatch(bigutilsr::tukey_mc_up(mah_dist),     error = function(e) NA_real_)
  )
  cand_MAH <- cand_MAH[is.finite(cand_MAH) & !is.na(cand_MAH)]
  
  cand_LOF <- c(
    tryCatch(bigutilsr::hist_out(lof_dist)$lim[2], error = function(e) NA_real_),
    tryCatch(bigutilsr::tukey_mc_up(lof_dist),     error = function(e) NA_real_)
  )
  cand_LOF <- cand_LOF[is.finite(cand_LOF) & !is.na(cand_LOF)]
  
  if (length(cand_MAH) == 0L && length(cand_LOF) == 0L) return(NULL)
  
  out_mah <- if (length(cand_MAH) > 0L) (mah_dist > min(cand_MAH)) else rep(FALSE, nrow(m))
  out_lof <- if (length(cand_LOF) > 0L) (lof_dist > min(cand_LOF)) else rep(FALSE, nrow(m))
  
  # Původní plná citlivost (OR)
  is_outlier <- out_mah | out_lof
  is_outlier[is.na(is_outlier)] <- FALSE
  
  return(d_clean$SEGMENT_ID[which(is_outlier)])
}

# Spuštění detekce
outlier_ids <- map(valid_biotops, detect_outliers) %>% unlist() %>% unique()
outlier_ids <- outlier_ids[!is.na(outlier_ids)]

cat("Identifikováno celkem", length(outlier_ids), "odlehlých segmentů.\n")

# 5. Aktualizace datasetu
df_model <- df_model %>% 
  mutate(outlier = SEGMENT_ID %in% outlier_ids)

df_outliers <- df_model %>% filter(outlier == TRUE) %>% select(SEGMENT_ID)

# 6. Prostorové připojení a export
cat("Načítám prostorová data z", file_shp, "...\n")
vmb <- st_read(file_shp, quiet = TRUE) %>% 
  mutate(SEGMENT_ID = as.character(SEGMENT_ID))

vmb_outliers <- vmb %>% inner_join(df_outliers, by = "SEGMENT_ID")

cat("Ukládám výsledný shapefile do:", out_shp, "\n")
st_write(vmb_outliers, out_shp, delete_layer = TRUE, quiet = TRUE)