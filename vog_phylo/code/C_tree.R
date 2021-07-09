# setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

#import tree
iqt <- read.newick(here("vog_phylo/data/align-trim-tree/iqtree1/sigmas_MafftEinsi.trim.iqtree"))

# list label data
d.iqt <- as_tibble(iqt)

#load alignment data
d.aln <- read_csv(here("data/sigmas_to_align.csv"))

# add data collected for bacteria from fearure tables
d.sp<- read_csv( here("data/bacterial_features.csv")) 
d.sp <- d.sp %>% 
  select(protein = product_accession, sp=sp, description=name, symbol)

d.aln <- rows_update(mutate(d.aln,symbol=""), d.sp, by = c("protein"))



# add meta data to tree tibble
d.iqt <- 
  left_join(d.iqt, d.aln, by = c("label" = "new.header"),)

# make informative label
d.iqt <- d.iqt %>% 
  mutate(tip.label = case_when(str_detect(description, regex("ecf", ignore_case = T)) ~ "ECF",
                               str_detect(sp, regex("subtilis", ignore_case = T)) ~ symbol,
                               TRUE ~ ""))
tree <- as.treedata(d.iqt)


# # unrooted trees
p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color=group), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)

ggsave(here("plots","sigma_all_unrooted.pdf"),p, height=10, width = 10)


# identify splitting node
p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color=group), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)+
  geom_text(aes(x=branch, label=node))  

ggsave(here("plots","sigma_nodeNUMS_unrooted.pdf"),p, height=10, width = 10)


# split tree at node 518
x <- groupClade(iqt,.node=518)
ggtree(x)+
  geom_tippoint(aes(color=group), size=1, shape=20)


# add split groups to main tree
d.x <- as_tibble(x) %>% 
  select(node, split = group)

d.iqt <- full_join(d.iqt, d.x)
tree <- as.treedata(d.iqt)

p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color = split), size=5, shape=21)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)


ggsave(here("plots","sigma_split_unrooted.pdf"),p, height=10, width = 10)

#need to keep split == 0

#split by group, and add ECF of B. subtilis, for reference.
d.keep <- d.iqt %>% 
  filter(split == 0 | str_detect(sp, regex("subtilis", ignore_case = T)))

keep <- d.keep %>%
  filter(!is.na(protein)) %>% 
  pull(protein)

##### Save results #####
if (!dir.exists(here("data", "reduced_set_to_align"))){
  dir.create(here("data", "reduced_set_to_align"))
}
# filter multifasta of hsearch results for filtering
library(seqinr)
sigma_fa <- read.fasta(here("data/sigmas_to_align.faa"),
                       seqtype =  "AA", whole.header = TRUE)

header.protein <- names(sigma_fa) %>% 
  str_remove("_bacteria") %>% 
  str_remove("_phage")

keep_fa <- sigma_fa[c(which(header.protein %in% keep))]

write.fasta(sequences = getSequence(keep_fa),
            names = getName(keep_fa),
            file.out = here("data", "reduced_set_to_align","sigmas_to_align.faa"))

# make filtered versions of meta data tables
d.faa <- read_csv(here("data/sigmas_to_align.csv"))
d.faa <- d.faa %>% 
  filter(protein %in% keep)
write_csv(d.faa, here("data", "reduced_set_to_align","sigmas_to_align.csv"))

d.sp<- read_csv( here("data/bacterial_features.csv"))
d.sp <- d.sp %>% 
  filter(product_accession %in% keep)
write_csv(d.sp, here("data", "reduced_set_to_align","bacterial_features.csv"))

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



as.treedata(d.iqt) %>% 
  ggtree( layout = 'equal_angle', aes(color = clade))

as.treedata(d.iqt) %>% 
  ggtree( layout = 'radial', aes(color = clade))+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = 1)

as.treedata(d.iqt) %>% 
  ggtree(aes(color = clade))+
  geom_tippoint(aes(color=group), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)
