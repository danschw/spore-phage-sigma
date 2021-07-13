#setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

#import tree
iqt <- read.iqtree(here("vog_phylo","data/curated_set_to_align/iqtree-fullALN-support/sigmas_MafftEinsi.aln.treefile"))

# iqt <- read.iqtree(here("vog_phylo","data/curated_set_to_align/iqtree1-support/sigmas_MafftEinsi.trim.contree"))
# iqt <- read.iqtree(here("vog_phylo","data/curated_set_to_align/iqtree1-support/sigmas_MafftEinsi.trim.treefile"))
# iqt <- read.newick(here("vog_phylo","data/align-trim-tree/sigmas_MafftEinsi.trim.treefile"))

# list label data
d.iqt <- as_tibble(iqt) %>% 
  mutate(group = case_when( str_detect(label, "bacteria") ~ "bacteria",
                            str_detect(label, "phage") ~ "phage",
                            TRUE ~ "NA")) %>% 
  mutate(protein.id = str_remove(label, "-.*"))


####load metadata for sequences ####
get_locus_tag <- function(protein_id){
  system(paste("wsl efetch -db protein -id", protein_id, "-format gb"), intern = T) %>% 
    tibble() %>% 
    filter(str_detect(., "locus_tag")) %>%
    filter(!str_detect(., "old_locus_tag")) %>%
    str_extract(., '(?<=").*?(?=")') %>% 
    return()
}

# load viral sigmas data from vog HMM analysis
# load(here("vogdb/data/vog_sigma_clean_Whost.RData"))# d.faa <- read_csv(here("data/sigmas_to_align.csv"))
# d.phage <- d.faa %>%
#   filter(protein %in% d.iqt$protein.id) %>% 
#   rename(protein.id = protein) %>% 
#   mutate(locus_tag = map_chr(protein.id, get_locus_tag))
# rm(d.faa)
# 
# # add 2 missing locus tags
# d.phage <- d.phage %>% 
#   mutate(locus_tag = case_when(protein.id == "YP_009031481.1" ~ "Bcp1_199",
#                                protein.id == "YP_009196983.1" ~ "CPT_Stills98",
#                                TRUE ~ locus_tag))
# 
# save(d.phage, file = here("vog_phylo/data/vog_sigma_curated.RData"))
load(file = here("vog_phylo/data/vog_sigma_curated.RData")) 

d.phage <- d.phage%>% 
  mutate(symbol = str_remove(locus_tag, ".*_"))%>% 
  mutate(lab_sp = str_remove(sp, ".* ")) %>% 
  mutate(lab_sp = str_remove(lab_sp, "vB_BsuM-")) %>% # for Goe3
  #remove phage name from locus tag
  mutate(symbol = str_remove(symbol, lab_sp)) %>% 
  mutate(symbol = if_else(lab_sp=="SPO1", description, symbol)) %>% 
  mutate(lab = paste(lab_sp, symbol, sep = "_"))



# load data collected for bacteria from fearure tables
d.bact<- read_csv( here("vog_phylo/data/bacterial_features.csv"))
d.bact <- d.bact %>%
  filter(product_accession %in% d.iqt$protein.id) %>% 
  rename(protein.id = product_accession, description=name) 

d.bact <- d.bact %>%
  mutate(lab_sp = case_when(str_detect(sp,"ubtilis") ~ "Bs",
                            str_detect(sp,"Clostridioides") ~ "Cd",
                            str_detect(sp,"Escherichia") ~ "Ec")) %>% 
  mutate(lab = paste(lab_sp, symbol, sep = "_"))

# add meta data to tree tibble
d.meta <-
  d.phage %>% 
    select( protein.id, description, sp, symbol, lab) %>% 
  bind_rows(., d.bact %>% select( protein.id, description, sp, symbol, lab))
  
d.iqt <- left_join(d.iqt, d.meta, by = "protein.id") 
  

#### root at base of ECF ####

# find node
d.iqt %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) %>% 
  as.treedata() %>% 
  ggtree()+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label = bs.label), color="blue", size=3, offset = .1)+
  geom_text(aes(x=branch, label=node), color="red", size=3) +
  ggsave(here("vog_phylo","plots","curated-sigma_nodeNUMS_unrooted.pdf"), height=10, width = 10)


# new root node
root_ecf <- 100

# assign root and add data
iqt <- root(iqt, node = root_ecf)

d.iqt <-  
  as.tibble(iqt) %>%
  mutate(group = if_else( str_detect(label, "bacteria"), "bacteria", "phage")) %>% 
  mutate(protein.id = str_remove(label, "-.*")) %>% 
  left_join(., d.meta, by = "protein.id")

# replot to verify re-rooting
d.iqt %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) %>% 
  as.treedata() %>% 
  ggtree()+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label = lab), color="blue", size=3, offset = .1)


#################
# function to mark phage only monophyletic clades
# go over all internal nodes, divide to clade and check if one side is phage only
# if yes mark as phage lade node

internal <- d.iqt %>% 
  filter(node > as.phylo(.) %>% Ntip()) %>% 
  pull(node)

#assign phage tip nodes to clade, and all other as empty vector
d.iqt <- d.iqt %>% 
  mutate(clade = case_when( str_detect(label, "bacteria") ~ "bacteria",
                            str_detect(label, "phage") ~ "phage",
                            TRUE ~ "bacteria"))

for (i in internal){
  x <- groupClade(iqt,.node=i)
  
  # add groups to main tree
  d.x <- as_tibble(x) %>% 
    select(node, cur.split = group) %>% 
    left_join(d.iqt,., by = "node")
  
  # get groups as charcater vectors
  g1 <- d.x %>% 
    filter(cur.split == 1) %>% 
    filter(! node %in% internal) %>% 
    pull(group)
  
  g2 <- d.x %>% 
    filter(cur.split == 0) %>% 
    filter(! node %in% internal) %>% 
    pull(group)
  
  # test if any of the groups is phage only
  if (length(g1)==0) next # avoid empty charcter returning TRUE
  if (length(g2)==0) next
  if (all(g1=="phage") | all(g2=="phage")){
    d.iqt$clade[d.iqt$node==i] <- "phage" # assign to phage only clade
  }
}


#################



# # first tree
ggtree(as.treedata(d.iqt) , aes(color = clade))+
  geom_tippoint(aes(shape=group), size=2, shape=20)+
  geom_tiplab(aes(label=lab), color="blue", size=3, offset = .1)

# # unrooted trees
p <-
  d.iqt %>% 
  mutate(support = case_when( #is.na(UFboot) ~ "NA",
    UFboot>= 90 ~ ">95%",
    UFboot>= 80 ~ ">80%",
    UFboot>= 70 ~ ">70%",
    TRUE ~ "<=70%")) %>% 
  as.treedata(.) %>% 
  ggtree(layout = 'equal_angle')+
  geom_nodepoint(aes(fill=support), size=1.5, shape=21)+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label=lab), color="blue", size=3, offset = .1)

ggsave(here("vog_phylo","plots","sigma_curated_unrooted.pdf"),p, height=10, width = 10)


as.treedata(d.iqt) %>% 
  ggtree( layout = 'equal_angle', aes(color = clade))

as.treedata(d.iqt) %>% 
  ggtree( layout = 'radial', aes(color = clade))+
  geom_tiplab(aes(label=lab), color="blue", size=3, offset = 1)

p <-as.treedata(d.iqt) %>% 
  ggtree(aes(color = clade))+
  # geom_tippoint(aes(fill=group), size=1, shape=20)+
  geom_tiplab(aes(label=lab), color="blue", size=3, offset = .1)+
  scale_color_manual(values = c("grey30", "blue"))

# ggsave(here("vog_phylo","plots","sigma_curated_rooted.pdf"),p, height=10, width = 10)


d.iqt %>% 
  mutate(support = case_when( #is.na(UFboot) ~ "NA",
                              UFboot>= 90 ~ ">95%",
                              UFboot>= 85 ~ ">85%",
                              UFboot>= 70 ~ ">70%",
                              TRUE ~ "")) %>% 
  as.treedata(.) %>% 
  ggtree(aes(color = clade))+#, layout = "fan")+
  geom_nodepoint(aes(fill=support), color=rgb(0,0,0,0), size=3, shape=22)+
  # geom_tippoint(aes(fill=group), size=1, shape=20)+
  geom_tiplab(aes(label=lab, color=clade), size=3, offset = .01, show.legend=F)+ 
  ggplot2::xlim(0, 6)+
  scale_color_manual(values = c("grey30", "blue"))+
  scale_fill_manual(values = c(rgb(0,0,0,0), grey.colors(3, rev = T) ))+
  theme(legend.position = "left")+
  ggsave(here("vog_phylo","plots","sigma_curated_rooted.pdf"), height=10, width = 8)



