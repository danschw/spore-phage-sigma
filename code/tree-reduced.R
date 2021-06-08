#setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

#import tree
iqt <- read.newick(here("data/reduced_set_to_align/multi-run-iqtree/sigmas_MafftEinsi.trim.treefile"))
# iqt <- read.newick(here("data/align-trim-tree/sigmas_MafftEinsi.trim.treefile"))

# list label data
d.iqt <- as_tibble(iqt)

#load alignment data
d.aln <- read_csv(here("data/reduced_set_to_align/sigmas_to_align.csv"))

# add data collected for bacteria from fearure tables
d.sp<- read_csv( here("data/reduced_set_to_align/bacterial_features.csv")) 
d.sp <- d.sp %>% 
  select(protein = product_accession, sp=sp, description=name, symbol)

d.aln <- rows_update(mutate(d.aln,symbol=""), d.sp, by = c("protein"))

# add meta data to tree tibble
d.iqt <- 
  left_join(d.iqt, d.aln, by = c("label" = "new.header"),)

#### root at base of ECF ####

# find node
d.iqt %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) %>% 
  as.treedata() %>% 
  ggtree()+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label = bs.label), color="blue", size=3, offset = .1)+
  geom_text(aes(x=branch, label=node)) 

# new root node
root_ecf <- 410

# assign root and add data
iqt <- root(iqt, node = root_ecf)

d.iqt <-  
  as.tibble(iqt) %>% 
  left_join(., d.aln, by = c("label" = "new.header"),)

# replot to verify re-rooting
d.iqt %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,"")) %>% 
  as.treedata() %>% 
  ggtree()+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label = bs.label), color="blue", size=3, offset = .1)

#### Plotting ####
# make informative label
d.iqt <- d.iqt %>% 
  mutate(tip.label = case_when(group == "phage" ~ sp,
                               group == "bacteria" ~ paste(sp, symbol, sep = "_"))) %>% 
  mutate(bs.label = if_else((group=="bacteria")&(str_detect(sp,"ubtilis")),
                            symbol,""))

#################
# function to mark phage only monophyletic clades
# go over all internal nodes, divide to clade and check if one side is phage only
# if yes mark as phage lade node

internal <- d.iqt %>% 
  filter(is.na(label)) %>% 
  pull(node)

#assign phage tip nodes to clade, and all other as empty vector
d.iqt <- d.iqt %>% 
  mutate(clade = if_else(group=="phage", "phage", "", missing = ""))

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
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = .1)

# # unrooted trees
p <-
  ggtree(as.treedata(d.iqt), layout = 'equal_angle')+
  geom_tippoint(aes(color=group), size=2, shape=20)+
  geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = .1)

ggsave(here("plots","sigma_reduced_unrooted.pdf"),p, height=10, width = 10)


as.treedata(d.iqt) %>% 
  ggtree( layout = 'equal_angle', aes(color = clade))

as.treedata(d.iqt) %>% 
  ggtree( layout = 'radial', aes(color = clade))+
  geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = 1)

p <-as.treedata(d.iqt) %>% 
  ggtree(aes(color = clade))+
  # geom_tippoint(aes(fill=group), size=1, shape=20)+
  geom_tiplab(aes(label=bs.label), color="blue", size=3, offset = .1)+
  scale_color_manual(values = c("grey30", "blue"))

ggsave(here("plots","sigma_reduced_rooted.pdf"),p, height=10, width = 10)
