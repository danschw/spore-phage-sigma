#setwd("/N/u/danschw/Carbonate/GitHub/spore-phage-sigma")
library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

#import tree
iqt <- read.newick(here("data/align-trim-tree/multi-iq/sigmas_MafftEinsi.trim.treefile"))

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
  # geom_point2(aes(subset=!isTip ,
  #                 fill=uf.boot),
  #             # fill=cut(uf.boot, c(0, 70, 100))),
  #             shape=21, size=1.5, color="white",stroke = 0) +
  geom_tippoint(aes(color=group), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)
# +#align = T, linetype = NULL)+
#   scale_fill_manual(values=c( "red"), guide='legend')+
#   scale_color_manual(values=c("grey","blue"))

ggsave(here("plots","sigma_all_unrooted.pdf"),p, height=10, width = 10)


# identify splitting node
p <-
  ggtree(tree, layout = 'equal_angle')+
  # geom_point2(aes(subset=!isTip ,
  #                 fill=uf.boot),
  #             # fill=cut(uf.boot, c(0, 70, 100))),
  #             shape=21, size=1.5, color="white",stroke = 0) +
  geom_tippoint(aes(color=group), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)+
  geom_text(aes(x=branch, label=node))  
# +#align = T, linetype = NULL)+
#   scale_fill_manual(values=c( "red"), guide='legend')+
#   scale_color_manual(values=c("grey","blue"))

ggsave(here("plots","sigma_nodeNUMS_unrooted.pdf"),p, height=10, width = 10)


# split tree at node 532
x <- groupClade(iqt,.node=532)
ggtree(x)+
  geom_tippoint(aes(color=group), size=1, shape=20)


# add groups to main tree
d.x <- as_tibble(x) %>% 
  select(node, split = group)

d.iqt <- full_join(d.iqt, d.x)
tree <- as.treedata(d.iqt)

p <-
  ggtree(tree, layout = 'equal_angle')+
  geom_tippoint(aes(color = split), size=1, shape=20)+
  geom_tiplab(aes(label=tip.label), color="blue", size=3, offset = .1)


ggsave(here("plots","sigma_split_unrooted.pdf"),p, height=10, width = 10)

#need to keep split == 0

#split by group, and add ECF of B. subtilis, for reference.
d.keep <- d.iqt %>% 
  filter(split == 0 | str_detect(sp, regex("subtilis", ignore_case = T)))

keep <- d.keep %>%
  filter(!is.na(protein)) %>% 
  pull(protein)

# # #####################################
# # # find clades of phage only proteins
# # phage.direct.parents <- 
# #   d.iqt%>%
# #   filter(group=="phage")%>%
# #   select(parent)%>%
# #   filter(!duplicated(parent))%>%
# #   as_vector()
# # 
# # d.iqt$clade <- NA
# # # assign phage node to phage clade
# # d.iqt <- 
# #   d.iqt%>%
# #   mutate(clade=if_else(group=="phage", group, "not"))#%>%
# #   # mutate(clade=na_if(clade,"not"))
# # 
# # # repeat agrregation 100 times unless done
# # for (r in 1:100){
# #   phage.parents <- c()
# #   for (i in phage.direct.parents){
# #     
# #     #filter root? parent and node are the same
# #     if (d.iqt$parent[d.iqt$node==i]==i) next
# #     
# #     #get direct children of current parent
# #     children <- child(iqt,i)
# #     
# #     #check if children nodes are all phage
# #     child.group <- 
# #       d.iqt%>%
# #       filter(node %in% children)%>%
# #       select(clade)%>%
# #       # +#align = T, linetype = NULL)+
#   scale_fill_manual(values=c( "red"), guide='legend')+
#   scale_color_manual(values=c("grey","blue"))

# as_vector()
# #     
# #     if(all(child.group=="phage", na.rm =T))
# #       phage.parents <- c(phage.parents,i)
# #   }
# #   #break when no more new upstream nodes are found
# #   if(is.null(phage.parents)) break
# #   #assign phage clade nodes
# #   d.iqt$clade[phage.parents] <- "phage"
# #   # find clades of phage only proteins
# #   phage.direct.parents <- 
# #     d.iqt[phage.parents,]%>%
# #     select(parent)%>%
# #     filter(!duplicated(parent))%>%
# #     as_vector()
# #   print(r)
# # }
# # 
# # d.iqt$clade[d.iqt$clade!="phage"] <- " "
# 
# # # # add bootstrap indicator
# #   d.iqt <-
# #     d.iqt%>%
# #     # mutate(uf.boot=as.numeric(label))%>%
# #   mutate(uf.boot=UFboot>70)
# #   d.iqt$uf.boot[!d.iqt$uf.boot] <- NA
# 
# # #add tip label
# #   tips <- length(iqt@phylo$tip.label)
# #   d.iqt$tip.lab <- NA
# #   d.iqt$tip.lab[d.iqt$node<=tips] <- d.iqt$label[d.iqt$node<=tips]
# 
# 
#   
# #root tree at sigN parent 
# # tree <- root(iqt,node=545)%>%
# #     as.treedata()
#   
# # tree <- 
# #   d.iqt%>%
# #   select(node,group,sp,protein,clade)%>%
# #   full_join(iqt, ., by='node')
# 
# tree <- as.treedata(d.iqt)
# 
# 
# # # unrooted trees
# p <-
#   ggtree(tree, layout = 'equal_angle')+
#   # geom_point2(aes(subset=!isTip ,
#   #                 fill=uf.boot),
#   #             # fill=cut(uf.boot, c(0, 70, 100))),
#   #             shape=21, size=1.5, color="white",stroke = 0) +
#   geom_tippoint(aes(color=group), size=1, shape=20)+
#   geom_tiplab(aes(label=symbol, color=group), size=3, offset = .1)
# # +#align = T, linetype = NULL)+
# #   scale_fill_manual(values=c( "red"), guide='legend')+
# #   scale_color_manual(values=c("grey","blue"))
# 
# p
# 
# 
# 
# 
# 
# 
# p <- 
# ggtree(tree, aes(color=clade),size=.5)+
#   # geom_point2(aes(subset=!isTip , 
#   #                 fill=uf.boot),
#   #             # fill=cut(uf.boot, c(0, 70, 100))), 
#   #             shape=21, size=1.5, color="white",stroke = 0) +
#   # geom_tippoint(aes(color=source), size=1, shape=20)+
#   geom_tiplab(aes(label=protein, color=clade), size=3, offset = .1)+#align = T, linetype = NULL)+
#     scale_fill_manual(values=c( "red"), guide='legend')+
#   scale_color_manual(values=c("pink","blue", "red"))+
#   theme(legend.position = "bottom")
#   
# p
# 
#   # ggsave(here("figures","ssp_tree.pdf"),p, height=10, width = 12)
#   
#  
# #   # tree to try out node rotations
# # p <- 
# # ggtree(tree)+
# #   geom_tippoint(aes(color=source), size=0.8)+
# #   geom_tiplab(aes(label=sp, color=source), size=3, offset = 2)+#align = T, linetype = NULL)+
# # geom_nodelab(aes(label=node))+
# #   scale_fill_manual(values=c("white", "blue"), guide='legend') 
# # p  
# # os <- offspring(tree,53)[1:2]
# # p <- p%>%flip(os[1],os[2])#%>%rotate(408)%>%rotate(409)%>%rotate(410)
# # # p+geom_nodelab(aes(label=node))
# # p 
# 
# # #####################
# # # unrooted trees
# p <-
#   ggtree(tree, aes(color=clade), layout = 'equal_angle')+
#   # geom_point2(aes(subset=!isTip ,
#   #                 fill=uf.boot),
#   #             # fill=cut(uf.boot, c(0, 70, 100))),
#   #             shape=21, size=1.5, color="white",stroke = 0) +
#   # geom_tippoint(aes(color=source), size=1, shape=20)+
#   geom_tiplab(aes(label=protein, color=group), size=3, offset = .1)+#align = T, linetype = NULL)+
#   scale_fill_manual(values=c( "red"), guide='legend')+
#   scale_color_manual(values=c("grey","blue"))
# 
# p
# # p <- p%>%rotate(405)%>%rotate(408)%>%rotate(409)%>%rotate(410)
# ggsave(here("figures","unrooted_ssp_tree.pdf"),p, height=10, width = 7)
# # 
# # 
# # #####################
# # # circular trees
# # # root.t <- rootnode(tree)
# # 
# p <-
#   ggtree(tree, aes(color=clade),size=.5, layout = 'circular',branch.length = "none")+
#   geom_point2(aes(subset=!isTip ,
#                   fill=uf.boot),
#               # fill=cut(uf.boot, c(0, 70, 100))),
#               shape=21, size=1.5, color="white",stroke = 0) +
#   # geom_tippoint(aes(color=source), size=1, shape=20)+
#   geom_tiplab(aes(label=tip.lab, color=clade), size=3, offset = 1)+#align = T, linetype = NULL)+
#   scale_fill_manual(values=c( "red"), guide='legend')+
#   scale_color_manual(values=c("grey","blue"))
#   ggsave(here("figures","ssp_circular_tree.pdf"),p, height=13, width = 13)
# # 
# # 
# #   ggtree(tree, aes(color=clade),size=.5, layout = 'circular',branch.length = "none")+
# #   geom_tippoint(aes(color=source), size=0.1)+
# #   geom_tiplab( aes(label=tip.lab,color=source),size=.2)+
# #   ggsave(here("circular_tree_labels.pdf"))
