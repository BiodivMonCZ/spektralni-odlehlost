# ==============================================================================
# compute_stats.R: (mean, median, sd)
# ==============================================================================

library(sf)
library(terra)
library(dplyr)
library(progress)

# Nastavení relativních cest vůči kořenu projektu
dir_interim <- "interim"
vmb_path    <- file.path(dir_interim, "VMB_filt.shp")
rast_path   <- file.path(dir_interim, "rast_stack.tif")

out_csv     <- file.path(dir_interim, "VMB_stats.csv")
out_gpkg    <- file.path(dir_interim, "VMB_stats.gpkg")

if (!file.exists(vmb_path)) stop("Nenalezen vstupní vektor: ", vmb_path)
if (!file.exists(rast_path)) stop("Nenalezen rastrový stack: ", rast_path)

# Načtení dat a negativní vnitřní buffer (-15 m)
cat("Načítám vektor:", vmb_path, "...\n")
vmb <- st_read(vmb_path, quiet = TRUE)
cat("Aplikuji vnitřní buffer -15 m na polygony...\n")
vmb <- st_buffer(vmb, -15)
vmb <- vmb[!st_is_empty(vmb), ]
cat("Zbývá platných polygonů po oříznutí bufferem:", nrow(vmb), "\n")

# Načtení rastrového stacku
cat("Načítám rastrový stack:", rast_path, "...\n")
r_stack <- rast(rast_path)

# Kontrola shody CRS
vmb_epsg <- st_crs(vmb)$epsg
rast_epsg <- terra::crs(r_stack, describe = TRUE)$code

if (is.na(vmb_epsg) || is.na(rast_epsg)) {
  warning("Jeden z CRS nemá definovaný EPSG kód. Kontrola není spolehlivá.")
} else if (vmb_epsg != rast_epsg) {
  warning(paste0("CRS se neshodují (VMB EPSG: ", vmb_epsg, " vs rastr EPSG: ", rast_epsg, ")."))
} else {
  cat(" CRS je shodný (oba EPSG:", vmb_epsg, ").\n")
}

# Ošetření duplicitních názvů vrstev v rastru
if (any(duplicated(names(r_stack)))) {
  names(r_stack) <- make.unique(names(r_stack), sep = "_")
}

# Iterace po REGION_ID (nebo po blocích)
region_ids <- unique(vmb$REGION_ID)
n_regions  <- length(region_ids)

cat("Zpracovávám", n_regions, "okrsků (REGION_ID)...\n")

stats_list <- vector("list", n_regions)

for (i in seq_along(region_ids)) {
  region_i <- region_ids[i]
  
  # Výpis stavu
  if (i == 1 || i %% 50 == 0 || i == n_regions) {
    pct <- round((i / n_regions) * 100, 1)
    cat(sprintf("[%s] Okrsek %d z %d (%.1f %%) - REGION_ID: %s\n", 
                format(Sys.time(), "%H:%M:%S"), i, n_regions, pct, region_i))
  }
  
  vmb_block <- vmb[vmb$REGION_ID == region_i, ]
  
  vmb_vect_block <- vect(vmb_block)
  
  # Výřez rastru podle rozsahu okrsku pro zrychlení výpočtu
  r_crop <- crop(r_stack, ext(vmb_vect_block))
  
  # Extrakce statistik
  stats <- terra::extract(
    r_crop,
    vmb_vect_block,
    fun = c("mean", "median", "sd"),
    na.rm = TRUE,
    touches = FALSE,
    ID = TRUE
  )
  
  # Lookup tabulka pro propojení ID polygonu
  lookup <- data.frame(
    ID = seq_len(nrow(vmb_block)),
    SEGMENT_ID = vmb_block$SEGMENT_ID
  )
  
  stats <- left_join(stats, lookup, by = "ID")
  
  if (nrow(stats) != nrow(vmb_block)) {
    warning(paste0("Nesoulad počtu řádků pro REGION_ID: ", region_i))
  }
  
  stats_list[[i]] <- stats
}

# Sjednocení a filtrace NA záznamů
stats_all <- bind_rows(stats_list)
stats_all <- dplyr::select(stats_all, -ID)

n_pred <- nrow(stats_all)
stats_all <- stats_all[complete.cases(stats_all), ]
cat("\nOdebráno", n_pred - nrow(stats_all), "segmentů kvůli chybějícím statistikám (NA).\n")
cat("Výsledný počet platných segmentů:", nrow(stats_all), "\n")

# Sloučení s atributy a export do interim/
cat("Zapisuji CSV výstup do:", out_csv, "...\n")
vmb_df <- st_drop_geometry(vmb)
vmb_out_df <- inner_join(vmb_df, stats_all, by = "SEGMENT_ID")
write.csv(vmb_out_df, out_csv, row.names = FALSE)

cat("Zapisuji GPKG výstup do:", out_gpkg, "...\n")
vmb_joined_sf <- inner_join(vmb, stats_all, by = "SEGMENT_ID")
st_write(vmb_joined_sf, out_gpkg, delete_layer = TRUE, quiet = TRUE)