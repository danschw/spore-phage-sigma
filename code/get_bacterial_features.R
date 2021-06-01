library(here)
library(tidyverse)


####################################################
# Get bacterial faa files specified by Burton et al. 2019
d.burton <- read_csv(here("data","Burton_S6.csv"))

#fix names to remove spaces
d.burton <- d.burton%>%
  mutate(organism_name=str_replace_all(organism_name," ","_"))

if (!dir.exists(here("data","bacteria_features/"))){
  dir.create(here("data","bacteria_features/"))
} 
setwd(here("data","bacteria_features/"))

d <- tibble()

for (i in 1:nrow(d.burton)){

  file.name <- paste0(str_remove(d.burton$ftp_path[i],".*/"),"_feature_table.txt.gz")
  
  ftp <- paste0(d.burton$ftp_path[i],"/",file.name)
  
  new.file.name <- paste0(d.burton$organism_name[i],"._features.tsv.gz")
  #download.file(ftp, destfile = new.file.name)
  
  d <- read_tsv(new.file.name) %>% 
    filter(`# feature` == "CDS") %>% 
    select(assembly, genomic_accession, product_accession, name, symbol) %>% 
    mutate(sp = d.burton$organism_name[i]) %>% 
    bind_rows(d,.)
}

#load alignment data
d.aln <- read_csv(here("data/sigmas_to_align.csv"))

# d.aln$group[d.aln$protein %in% d$product_accession] %>% table()
# d.aln$group[!d.aln$protein %in% d$product_accession] %>% table()

d.sp <- d %>% 
  filter(product_accession %in% d.aln$protein)

write_csv(d.sp, here("data/bacterial_features.csv"))
