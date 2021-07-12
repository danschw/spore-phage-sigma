#setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

#import tree
iqt <- read.iqtree(here("vog_phylo","data/reduced_set_to_align/iqtree1-support/sigmas_MafftEinsi.trim.treefile"))
# iqt <- read.newick(here("vog_phylo","data/align-trim-tree/sigmas_MafftEinsi.trim.treefile"))

# list label data
d.iqt <- as_tibble(iqt) %>% 
  mutate(group = case_when( str_detect(label, "bacteria") ~ "bacteria",
                            str_detect(label, "phage") ~ "phage",
                            TRUE ~ "NA")) %>% 
  mutate(protein.id = str_remove(label, "-.*"))


####load metadata for sequences ####

# load viral sigmas data from vog HMM analysis
load(here("vogdb/data/vog_sigma_clean_Whost.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
d.phage <- d.faa %>%
  filter(protein %in% d.iqt$protein.id) %>% 
  rename(protein.id = protein)
rm(d.faa)

# load data collected for bacteria from fearure tables
d.bact<- read_csv( here("vog_phylo/data/bacterial_features.csv"))
d.bact <- d.bact %>%
  filter(product_accession %in% d.iqt$protein.id) %>% 
  rename(protein.id = product_accession, description=name)

# add meta data to tree tibble
d.meta <-
  d.phage %>% 
  select( protein.id, description, sp) %>% 
  mutate(symbol=NA) %>% 
  bind_rows(., d.bact %>% select( protein.id, description, sp, symbol))

d.iqt <- left_join(d.iqt, d.meta, by = "protein.id")

#######

# siphos with 3 sigmas
siphos3 <- d.phage %>% 
  filter(str_detect(viral.family, "Sipho")) %>% 
  group_by(sp) %>% 
  summarise(n=n()) %>% 
  slice_max(n)


d.pick <- 
  d.iqt %>% 
  filter((group=="bacteria")&(str_detect(sp,"ubtilis")) |
           (group=="bacteria")&(str_detect(sp,"Clostridioides")) |
           (group=="bacteria")&(str_detect(sp,"Escherichia")) |
           (group=="phage")&(str_detect(sp,"SPO1")) | 
           (group=="phage")&(str_detect(sp,"SP-10"))|
           (group=="phage")&(str_detect(sp,"Goe3"))|
           (group=="phage")&(str_detect(sp,"Eldridge"))|
           (group=="phage")&(str_detect(sp,"Moonbeam"))|
           (group=="phage")&(str_detect(sp,"Fah"))|
           (group=="phage")&(str_detect(sp,"Bcp"))|
           sp %in% siphos3$sp)

         