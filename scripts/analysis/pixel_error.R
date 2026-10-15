# ==============================================================================
# pixel_error.R
# ==============================================================================

library(terra)
library(robust)      # pro tukey_mc_up
library(robustbase)  # pro covOGK
library(bigutilsr)   # pro hist_out
library(dplyr)

# 1. Nastavení cest
dir_interim <- "interim"
dir_outputs <- "outputs"

file_stack <- file.path(dir_interim, "rast_stack.tif")
file_vmb   <- file.path(dir_interim, "VMB_stats.gpkg")
out_tif    <- file.path(dir_outputs, "VMB_outliers_pixel.tif")

if (!file.exists(file_stack)) stop("Chybí vstupní rastrový stack: ", file_stack)
if (!file.exists(file_vmb)) stop("Chybí vstupní vektory statistik: ", file_vmb)

# 2. Načtení dat
cat("Načítám data...\n")
rast_stack <- terra::rast(file_stack)
vmb        <- terra::vect(file_vmb)

# Všechny stávající vrstvy ve stacku slouží jako prediktory
rast_vars <- names(rast_stack)

# 3. Nalezení platných biotopů (min. 31 polygonů) a číselné kódování
vmb_df <- as.data.frame(vmb)
biotop_counts <- table(vmb_df$BIOTOP_SEZ)
valid_biotopy <- names(biotop_counts)[biotop_counts >= 31]

cat("Nalezeno", length(valid_biotopy), "kategorií splňujících limit (>= 31 polygonů).\n\n")

# Explicitní ID pro párování rastru
biotop_lut <- data.frame(
  id = seq_along(valid_biotopy),
  biotop = valid_biotopy,
  stringsAsFactors = FALSE
)
vmb$biotop_id <- biotop_lut$id[match(vmb$BIOTOP_SEZ, biotop_lut$biotop)]

# 4. Trénink PCA a OGK pro každý platný biotop
models <- list()
n_valid <- length(valid_biotopy)

cat("--- Zahájení trénování modelů pro jednotlivé biotopy ---\n")

for (idx_b in seq_along(valid_biotopy)) {
  biotop <- valid_biotopy[idx_b]
  b_id   <- biotop_lut$id[idx_b]
  
  cat(sprintf("[%d/%d] Biotop '%s' (ID %d): ", idx_b, n_valid, biotop, b_id))
  
  d_poly <- vmb[vmb$biotop_id == b_id & !is.na(vmb$biotop_id), ]
  
  # Nalezení buněk mřížky pro daný biotop
  cell_res <- tryCatch(
    terra::cells(rast_stack[[1]], d_poly, touches = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(cell_res) || length(cell_res) == 0) {
    cat("PŘESKOČENO (nenalezeny žádné buňky na gridu).\n")
    next
  }
  
  # Získání unikátních indexů buněk
  cells_found <- if (is.matrix(cell_res)) {
    if ("cell" %in% colnames(cell_res)) cell_res[, "cell"] else cell_res[, ncol(cell_res)]
  } else {
    cell_res
  }
  
  cells_biotop <- unique(cells_found[!is.na(cells_found)])
  n_pix_total  <- length(cells_biotop)
  
  if (n_pix_total < 10) {
    cat(sprintf("PŘESKOČENO (nedostatek platných buněk: %d).\n", n_pix_total))
    next
  }
  
  # Reprezentativní výběr (max 10 000 buněk)
  set.seed(42)
  if (n_pix_total > 10000) {
    cells_sample <- sort(sample(cells_biotop, 10000))
  } else {
    cells_sample <- sort(cells_biotop)
  }
  
  # Čtení pouze vybraných buněk přímo z disku
  pix_data <- rast_stack[cells_sample]
  pix_clean <- na.omit(pix_data[, rast_vars, drop = FALSE])
  n_pix_used <- nrow(pix_clean)
  
  if (n_pix_used < 10) {
    cat("PŘESKOČENO (nedostatek kompletních spektrálních dat bez NA).\n")
    next
  }
  
  X <- pix_clean
  
  # Odstranění sloupců s nulovým rozptylem
  X_var <- apply(X, 2, var, na.rm = TRUE)
  X_model <- X[, X_var > 1e-6, drop = FALSE]
  
  if (ncol(X_model) < 2) {
    cat("PŘESKOČENO (nedostatek sloupců s rozptylem).\n")
    next
  }
  
  # PCA
  pca <- tryCatch(prcomp(X_model, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) {
    cat("PŘESKOČENO (selhalo PCA).\n")
    next
  }
  
  # Fixní výběr max. 10 komponent
  n_comp <- min(10, ncol(pca$x))
  m_pix <- pca$x[, 1:n_comp, drop = FALSE]
  
  # Robustní OGK
  ogk <- tryCatch(robustbase::covOGK(m_pix, sigmamu = robustbase::scaleTau2), error = function(e) NULL)
  if (is.null(ogk)) {
    cat("PŘESKOČENO (selhal výpočet covOGK).\n")
    next
  }
  
  rd_center <- as.numeric(ogk$center)
  rd_cov    <- ogk$cov
  
  # Odmocněná Mahalanobisova vzdálenost D
  D_pix <- sqrt(mahalanobis(m_pix, center = rd_center, cov = rd_cov))
  
  # Prahování
  thr_hist <- tryCatch(bigutilsr::hist_out(D_pix)$lim[2], error = function(e) NA_real_)
  thr_tuk  <- tryCatch(robust::tukey_mc_up(D_pix),        error = function(e) NA_real_)
  
  cand <- c(thr_hist, thr_tuk)
  cand <- cand[is.finite(cand) & !is.na(cand)]
  
  if (length(cand) == 0L) {
    cat("PŘESKOČENO (selhal výpočet prahu odlehlosti).\n")
    next
  }
  
  selected_thr <- min(cand)
  
  # Uložení modelu pod ID
  models[[as.character(b_id)]] <- list(
    pca       = pca,
    rd_center = rd_center,
    rd_cov    = rd_cov,
    k         = n_comp,
    mah_thr   = selected_thr,
    pix_vars  = names(X_model)
  )
  
  cat(sprintf("OK | Vzorkováno: %d/%d px | PC: %d | Práh D: %.3f\n", 
              n_pix_used, n_pix_total, n_comp, selected_thr))
}

cat("------------------------------------------------------\n")
if (length(models) == 0L) stop("Nebyl natrénován žádný model.")
cat(sprintf("Celkem úspěšně natrénováno modelů: %d z %d.\n\n", length(models), n_valid))

# 5. Rasterizace polygonů biotopů
cat("Rasterizuji polygony biotopů..\n")
biotop_rast <- terra::rasterize(vmb, rast_stack[[1]], field = "biotop_id", background = NA, touches = FALSE)
names(biotop_rast) <- "biotop_id"
rast_stack <- c(rast_stack, biotop_rast)

# 6. Bloková detekce odlehlostí
out_rast <- rast_stack[[1]]
values(out_rast) <- NA
terra::writeStart(out_rast, filename = out_tif, overwrite = TRUE, datatype = "INT1U", NAflag = 255)

bs <- terra::blocks(rast_stack)
cat("Spouštím blokovou predikci odlehlostí (", bs$n, "dlaždic )...\n")

failed_blocks <- integer(0)

for (i in seq_len(bs$n)) {
  if (i == 1 || i %% 10 == 0 || i == bs$n) {
    cat(sprintf("Zpracovávám dlaždici %d z %d...\n", i, bs$n))
  }
  
  success  <- FALSE
  attempts <- 0
  
  while (!success && attempts < 3) {
    attempts <- attempts + 1
    tryCatch({
      
      block <- terra::values(rast_stack, row = bs$row[i], nrows = bs$nrows[i], dataframe = TRUE)
      if (is.null(block) || nrow(block) == 0) stop("Prázdná dlaždice.")
      
      biotop_vec <- as.character(block[, "biotop_id"])
      outlier_block <- rep(NA_integer_, nrow(block))
      
      for (b_id_val in na.omit(unique(biotop_vec))) {
        b_id_chr <- as.character(b_id_val)
        if (!b_id_chr %in% names(models)) next
        
        model <- models[[b_id_chr]]
        
        # Výběr pixelů pro daný biotop a kompletní data
        idx <- which(biotop_vec == b_id_chr & complete.cases(block[, model$pix_vars, drop = FALSE]))
        if (length(idx) == 0) next
        
        d_pix <- block[idx, model$pix_vars, drop = FALSE]
        
        # Predikce PCA a odmocněná Mahalanobisova vzdálenost D
        pr <- tryCatch(predict(model$pca, newdata = d_pix), error = function(e) NULL)
        
        if (!is.null(pr) && ncol(pr) >= model$k) {
          pcs_pix <- pr[, 1:model$k, drop = FALSE]
          
          D_pix <- sqrt(mahalanobis(pcs_pix, center = model$rd_center, cov = model$rd_cov))
          outlier_block[idx] <- as.integer(D_pix > model$mah_thr)
        }
      }
      
      terra::writeValues(out_rast, outlier_block, bs$row[i], nrows = bs$nrows[i])
      success <- TRUE
      
    }, error = function(e) {
      if (attempts == 3) warning(paste("Dlaždice", i, "selhala:", conditionMessage(e)))
    })
    
    if (!success && attempts == 3) failed_blocks <- c(failed_blocks, i)
  }
}

out_rast <- terra::writeStop(out_rast)

if (length(failed_blocks) > 0) {
  cat("Selhaly dlaždice:", paste(failed_blocks, collapse = ", "), "\n")
} else {
  cat("Predikce byla úspěšně dokončena. Výstup uložen do:", out_tif, "\n")
}