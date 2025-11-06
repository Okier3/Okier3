# Load libraries
library(ggplot2)

# Read the data
setwd("D:/KE-TZ BOUNDARY PHASE 3/PROJECT [001]/EXCEL/R/Eg/proj_01")

ASTRO_file <- file.choose()
EGM96_file <- file.choose()
EGM2008_file <- file.choose()
EIGEN_file <- file.choose()
GECO_file <- file.choose()
SGG_file <- file.choose()

# Read CSV files
ASTRO_data <- read.csv(ASTRO_file, stringsAsFactors = FALSE)[, c("BPS_Name", "Xi_arcsec", "Eta_arcsec")]
EGM96_data <- read.csv(EGM96_file, stringsAsFactors = FALSE)[, c("BPS_Name", "Xi_arcsec", "Eta_arcsec")]
EGM2008_data <- read.csv(EGM2008_file, stringsAsFactors = FALSE)[, c("BPS_Name", "Xi_arcsec", "Eta_arcsec")]
EIGEN_data <- read.csv(EIGEN_file, stringsAsFactors = FALSE)[, c("BPS_Name", "Xi_arcsec", "Eta_arcsec")]
GECO_data <- read.csv(GECO_file, stringsAsFactors = FALSE)[, c("BPS_Name", "Xi_arcsec", "Eta_arcsec")]
SGG_data <- read.csv(SGG_file, stringsAsFactors = FALSE)[, c("BPS_Name", "Xi_arcsec", "Eta_arcsec")]

# Preview data
head(EGM96_data)
# column to identify source
ASTRO_data$Source   <- "ASTRO"
EGM96_data$Source   <- "EGM96"
EGM2008_data$Source <- "EGM2008"
EIGEN_data$Source   <- "EIGEN"
GECO_data$Source    <- "GECO"
SGG_data$Source     <- "SGG"
combined <- rbind(ASTRO_data, EGM96_data, EGM2008_data, EIGEN_data, GECO_data, SGG_data)
head(combined)

# Ensure PointName is treated as a factor so it stays in the right order
combined$BPS_Name <- factor(combined$BPS_Name, levels = unique(ASTRO_data$BPS_Name))

# Plot north south values
ggplot(combined, aes(x = BPS_Name, y = Eta_arcsec, color = Source, group = Source)) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_color_manual(values = c(
    "ASTRO"   = "red",
    "EGM96"   = "blue",
    "EGM2008" = "black",
    "EIGEN"   = "green",
    "GECO"    = "purple",
    "SGG"     = "orange"
  )) +
  # vertical partition lines for each point
  geom_vline(xintercept = seq_along(unique(combined$BPS_Name)), 
             color = "grey50", linetype = "dashed", size = 0.3) +
  # horizontal partition lines at Y breaks
  geom_hline(yintercept = pretty(combined$Eta_arcsec),
             color = "grey50", linetype = "dashed", size = 0.3) +
  labs(title = "Comparison of North-South component",
       x = "Point Name",
       y = "Eta (ξ)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", size = 0.8), # darker X & Y axis lines
        panel.grid.major.x = element_blank() # hide default grid, keep only partitions
        )

# Plot east West values
ggplot(combined, aes(x = BPS_Name, y = Xi_arcsec, color = Source, group = Source)) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_color_manual(values = c(
    "ASTRO"   = "red",
    "EGM96"   = "blue",
    "EGM2008" = "black",
    "EIGEN"   = "green",
    "GECO"    = "purple",
    "SGG"     = "orange"
  )) +
  # vertical partition lines for each point
  geom_vline(xintercept = seq_along(unique(combined$BPS_Name)), 
             color = "grey50", linetype = "dashed", size = 0.3) +
  # horizontal partition lines at Y breaks
  geom_hline(yintercept = pretty(combined$Eta_arcsec),
             color = "grey50", linetype = "dashed", size = 0.3) +
  labs(title = "Comparison of East-West component",
       x = "Point Name",
       y = "Xi (η)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", size = 0.8), # darker X & Y axis lines
        panel.grid.major.x = element_blank() # hide default grid, keep only partitions
        )
