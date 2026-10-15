# ==============================================================================
# vmb_filt.R
# ==============================================================================

library(sf)
library(dplyr)

# Definice cest
input_shp   <- "data/VMB/Aktualni_Segment.shp"
output_shp  <- "interim/VMB_filt.shp"

# Kontrola existence vstupního souboru
if (!file.exists(input_shp)) {
  stop("Vstupní soubor nebyl nalezen: ", input_shp)
}

# Načtení vrstvy
cat("Načítám:", input_shp, "...\n")
vmb <- st_read(input_shp, quiet = TRUE)
n_puvodni <- nrow(vmb)
cat("Načteno segmentů:", n_puvodni, "\n")

# Filtrace podle sloupce FSB
if ("FSB" %in% names(vmb)) {
  vmb <- vmb %>% 
    filter(!is.na(FSB)) %>% 
    filter(!(FSB %in% c("X", "moz.", "-")))
  
  cat("Odfiltrováno segmentů:", n_puvodni - nrow(vmb), "\n")
} else {
  warning("Chyba! Filtrace nebyla provedena.")
}

# Transformace souřadnicového systému do UTM 33N (EPSG:32633)
cat("Transformuji do EPSG:32633 \n")
vmb <- st_transform(vmb, crs = 32633)

# Uložení výsledku do složky interim
cat("Ukládám filtrovaný shapefile do:", output_shp, "...\n")
st_write(vmb, output_shp, delete_layer = TRUE, quiet = TRUE)
