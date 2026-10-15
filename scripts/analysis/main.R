# ==============================================================================
# Spektrální odlehlost
# ==============================================================================

library(callr)

project_root <- "C:/Users/jezek/Desktop/Data/spektralni_odlehlost_2026" # Nastav cestu
setwd(project_root)

# Definice cest
dir_scripts    <- file.path(project_root, "spektralni_odlehlost_R")
dir_data       <- file.path(project_root, "data")

# Inputs
dir_vmb <- file.path(dir_data, "VMB")
dir_s2  <- file.path(dir_data, "s2gm_copernicus")
dirs_s1 <- list.dirs(dir_data, recursive = FALSE, full.names = TRUE)
dirs_s1 <- dirs_s1[grepl("/S1_", dirs_s1)] # Složky dle prefixu S1_

# Meziprodukty a finální výstupy
dir_interim <- file.path(project_root, "interim")
dir_outputs <- file.path(project_root, "outputs")

# Založení složek
dir.create(dir_interim, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_outputs, recursive = TRUE, showWarnings = FALSE)

# Spouštěcí funkce
run_step <- function(script_name, step_title) {
  script_path <- file.path(dir_scripts, script_name)
  
  if (!file.exists(script_path)) {
    stop("Skript nebyl nalezen: ", script_path)
  }
  
  cat("\n======================================================================\n")
  cat(" SPUŠTĚNÍ:", step_title, "\n")
  cat(" Skript:  ", script_name, "\n")
  cat(" Start:   ", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("======================================================================\n\n")

  callr::rscript(
    script = script_path,
    wd = project_root,
    show = TRUE,
    fail_on_status = TRUE
  )
  
  cat("Dokončen:", step_title, "v", format(Sys.time(), "%H:%M:%S"), "\n")
}

# Spuštění kroků
run_step("vmb_filt.R", "1/5 Filtrace VMB")
run_step("sentinel_to_stack.R", "2/5 Mozaikování a tvorba rast_stack")
run_step("compute_stats.R", "3/5 Výpočet zonálních statistik")
run_step("polygon_error.R", "4/5 Polygonová detekce odlehlostí")
run_step("pixel_error.R", "5/5 Pixelová detekce odlehlostí")

cat("Hotovo")
