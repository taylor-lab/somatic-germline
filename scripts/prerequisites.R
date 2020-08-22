

# Required packages -----------------------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(plyr)
  library(ggsignif)
  library(binom)
  library(stringi)
  library(ggrepel)
  library(grid)
  library(readxl)
  library(here)
  library(survival)
  library(survminer)
  library(forestmodel)
})

# Load data -------------------------------------------------------------------------------------------------------
load('../data/supp_data_github.Rdata')

#######
# Set plotting theme ----------------------------------------------------------------------------------------------
theme_set(theme_bw())
theme_update(
  text = element_text(family = 'ArialMT', color = 'black', size = 12),
  axis.text = element_text(family = 'ArialMT', color = 'black', size = 12),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  legend.background = element_blank(),
  legend.key = element_blank()
)

