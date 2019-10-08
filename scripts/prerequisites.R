

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
})

# Load data -------------------------------------------------------------------------------------------------------
germline_cnv <- fread('../data/germline_cnv_V2.txt')
germline_variants <- fread('../data/germline_mutations_V3.txt')

cancer_types = read_xlsx('../data/supp_tables/Table S1 - Cancer Types in cohort.xlsx')
gene_variant_level_pen = read_xlsx('../data/supp_tables/Table S2 - Gene and variant-level penetrance assignments.xlsx') %>% janitor::clean_names()
gene_cancer_type_assoc = read_xlsx('../data/supp_tables/Table S5 - Gene and cancer type associations.xlsx') %>% janitor::clean_names()
clinical_age_dx = read_xlsx('../data/supp_tables/Table Sx - Clinical.xlsx')

load('../data/supp_data_github.Rdata')

germline_variants <-
  germline_variants %>%
  left_join(gene_variant_level_pen %>% select(gene, type) %>% unique, by=c("Hugo_Symbol" = "gene")) %>% 
  mutate(is_tumor_suppressor = ifelse(type == 'TS', T, F)) %>%
  mutate(sample_facets_qc_passed = ifelse(Normal_Sample %in% samples_facets_qc_passed, T, F)) %>%
  mutate(Hugo_Symbolp = paste0(Hugo_Symbol, ':', substr(penetrance,1,1)))

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

