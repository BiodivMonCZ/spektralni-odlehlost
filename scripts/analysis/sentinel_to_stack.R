# ==============================================================================
# sentinel_to_stack.R
# ==============================================================================

library(terra)

# 1. Nastavení cest
dir_data    <- "data"
dir_interim <- "interim"
s2_root     <- file.path(dir_data, "s2gm_copernicus")

if (!dir.exists(s2_root)) {
  stop("Vstupní soubor nebyl nalezen: ", s2_root)
}

# ČÁST 1: Sentinel-2 mozaikování z VRT
vrt_dirs <- unique(dirname(list.files(s2_root, pattern = "\\.vrt$", full.names = TRUE, recursive = TRUE)))

if (length(vrt_dirs) == 0) {
  stop("Nebyly nalezeny žádné .vrt soubory.")
}

term_dirs <- unique(dirname(vrt_dirs))

cat("Nalezeny termíny S2:\n")
for (td in term_dirs) cat(" -", td, "\n")
cat("\n")

find_tile_dirs <- function(term_dir) {
  cand <- list.dirs(term_dir, recursive = TRUE, full.names = TRUE)
  cand[grepl("^tile\\d+$|^tile$|tile\\d+|(^|[\\/])tile([^a-zA-Z]|$)", basename(cand), ignore.case = TRUE)]
}

for (term_dir in term_dirs) {
  term_name <- basename(term_dir)
  tile_dirs <- find_tile_dirs(term_dir)
  tile_dirs <- tile_dirs[dir.exists(tile_dirs)]
  
  s2_files <- character(0)
  if (length(tile_dirs) > 0) {
    s2_files <- unlist(lapply(tile_dirs, function(d) {
      list.files(path = d, pattern = "^B.*\\.tif$", full.names = TRUE, recursive = TRUE)
    }))
  }
  
  bands_from_tif <- character(0)
  if (length(s2_files) > 0) {
    bands_from_tif <- unique(sub("^.*(B0[2-8]|B1[1-2]|B8A)\\.tif$", "\\1", basename(s2_files)))
  }
  
  vrt_files <- list.files(term_dir, pattern = "^(B0[2-8]|B1[1-2]|B8A)\\.vrt$", full.names = TRUE, recursive = TRUE)
  bands_from_vrt <- unique(sub("\\.vrt$", "", basename(vrt_files)))
  
  bands <- unique(c(bands_from_tif, bands_from_vrt))
  bands <- bands[nzchar(bands)]
  
  if (length(bands) == 0) {
    cat("! Termín", term_name, ": nenalezen žádný band. Přeskakuji.\n")
    next
  }
  
  cat("Termín", term_name, ": zpracování pásem →", paste(sort(bands), collapse = ", "), "\n")
  
  for (band in sort(bands)) {
    band_vrt <- vrt_files[basename(vrt_files) == paste0(band, ".vrt")]
    if (length(band_vrt) == 0) next
    
    vrt_file <- band_vrt[1]
    mosaic_s2 <- try(rast(vrt_file), silent = TRUE)
    if (inherits(mosaic_s2, "try-error")) {
      cat(" ! Chyba načtení VRT:", vrt_file, "\n")
      next
    }
    
    output_file <- file.path(dir_interim, paste0("S2_", term_name, "_", band, "_mosaic.tif"))
    writeRaster(mosaic_s2, output_file, overwrite = TRUE)
    cat("S2 mozaika uložena:", basename(output_file), "\n")
  }
}

# ČÁST 2: Sentinel-1 mozaikování
mosaic_s1 <- function(folder_path, output_file) {
  tile_files <- list.files(path = folder_path, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  
  if (length(tile_files) == 0) {
    cat("Žádné TIF soubory v:", folder_path, "\n")
    return(NULL)
  }
  
  tiles <- lapply(tile_files, rast)
  mosaic_raster <- do.call(mosaic, c(tiles, fun = mean))
  
  output_path <- file.path(dir_interim, output_file)
  writeRaster(mosaic_raster, output_path, overwrite = TRUE)
  cat("S1 mozaika uložena:", basename(output_path), "\n")
}

s1_folders <- list.dirs(dir_data, full.names = TRUE, recursive = FALSE)
s1_folders <- s1_folders[grepl("[/\\\\]S1_", s1_folders)]

cat("\nNalezeno", length(s1_folders), "složek Sentinel-1 ke zpracování.\n")

for (folder in s1_folders) {
  folder_name <- basename(folder)
  output_file <- paste0(folder_name, "_mosaic.tif")
  cat(">>> Zpracovávám:", folder_name, "\n")
  mosaic_s1(folder_path = folder, output_file = output_file)
}

# ČÁST 3: Převzorkování a složení rast_stack.tif
cat("\n--- Tvorba společného stacku ---\n")

final_mosaics_paths <- list.files(
  path = dir_interim,
  pattern = "_mosaic\\.tif$",
  full.names = TRUE
)

# Automatické nalezení referenčního B02 pásma z vytvořených S2 mozaik
ref_candidates <- final_mosaics_paths[grepl("_B02_mosaic\\.tif$", final_mosaics_paths)]

if (length(ref_candidates) == 0) {
  stop("Nebyl nalezen žádný referenční raster B02")
}

ref_rast_path <- ref_candidates[1]
ref_rast <- rast(ref_rast_path)
cat("Použit referenční rastr pro převzorkování:", basename(ref_rast_path), "\n")

aligned_rasters <- list()

for (mosaic_file in final_mosaics_paths) {
  r <- rast(mosaic_file)
  
  # Pokud se geometrie/rozlišení liší, převzorkujeme na mřížku B02
  if (!compareGeom(r, ref_rast, stopOnError = FALSE)) {
    r <- resample(r, ref_rast, method = "bilinear")
    cat(" Převzorkováno:", basename(mosaic_file), "\n")
  } else {
    cat(" Shodná geometrie (bez převzorkování):", basename(mosaic_file), "\n")
  }
  
  aligned_rasters[[length(aligned_rasters) + 1]] <- r
}

if (length(aligned_rasters) > 0) {
  rast_stack <- rast(aligned_rasters)
  output_stack_file <- file.path(dir_interim, "rast_stack.tif")
  
  cat("\nZapisuji finální rast_stack.tif \n")
  writeRaster(rast_stack, output_stack_file, overwrite = TRUE)
  cat("Hotovo! Finální stack uložen do:", output_stack_file, "\n")
} else {
  stop("Chyba: Žádné mozaiky nebyly načteny pro vytvoření stacku.")
}