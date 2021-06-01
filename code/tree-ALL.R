library(here)
library(tidyverse)
library(ggtree) #https://yulab-smu.top/treedata-book/
library(treeio)

#import tree
iqt <- read.iqtree(here("data/align/iqtree/201221214318/out.trim.treefile"))

# list label data
d.iqt <- as_tibble(iqt)

#parse out name
d.iqt <- 
  d.iqt%>%
  mutate(lab=str_replace(label,"_",";"))%>%
  mutate(lab=str_replace(lab,"_",";"))%>%
  separate(lab, into = c("sp","locus","protein"),sep = ";")


# add phage label to tree based on phage hits
d.iqt <-
  d.iqt%>%
  mutate(phage =(str_detect( label,"phage")|str_detect(label,"virus")))%>%
  mutate(source=if_else(phage,"phage", "bacteria"))


#####################################
# find clades of phage only proteins
phage.direct.parents <- 
  d.iqt%>%
  filter(source=="phage")%>%
  select(parent)%>%
  filter(!duplicated(parent))%>%
  as_vector()

d.iqt$clade <- NA
# assign phage node to phage clade
d.iqt <- 
  d.iqt%>%
  mutate(clade=if_else(source=="phage", source, "not"))#%>%
  # mutate(clade=na_if(clade,"not"))

# repeat agrregation 100 times unless done
for (r in 1:100){
  phage.parents <- c()
  for (i in phage.direct.parents){
    #get direct children of current parent
    children <- child(iqt,i)
    
    #check if children nodes are all phage
    child.source <- 
      d.iqt%>%
      filter(node %in% children)%>%
      select(clade)%>%
      as_vector()
    
    if(all(child.source=="phage"))
      phage.parents <- c(phage.parents,i)
  }
  #break when no more new upstream ndes are found
  if(is.null(phage.parents)) break
  #assign phage clade nodes
  d.iqt$clade[phage.parents] <- "phage"
  # find clades of phage only proteins
  phage.direct.parents <- 
    d.iqt[phage.parents,]%>%
    select(parent)%>%
    filter(!duplicated(parent))%>%
    as_vector()
  print(r)
}

d.iqt$clade[d.iqt$clade!="phage"] <- " "





# # add bootstrap indicator
  d.iqt <-
    d.iqt%>%
    # mutate(uf.boot=as.numeric(label))%>%
  mutate(uf.boot=UFboot>70)
  d.iqt$uf.boot[!d.iqt$uf.boot] <- NA

#add tip label
  tips <- length(iqt@phylo$tip.label)
  d.iqt$tip.lab <- NA
  d.iqt$tip.lab[d.iqt$node<=tips] <- d.iqt$label[d.iqt$node<=tips]


  
#root tree at sigN parent 
# tree <- root(iqt,node=545)%>%
#     as.treedata()
  
tree <- 
  d.iqt%>%
  select(node,phage,source,sp,locus,protein,uf.boot,tip.lab,clade)%>%
  full_join(iqt, ., by='node')



p <- 
ggtree(tree, aes(color=clade),size=.5)+
  geom_point2(aes(subset=!isTip , 
                  fill=uf.boot),
              # fill=cut(uf.boot, c(0, 70, 100))), 
              shape=21, size=1.5, color="white",stroke = 0) +
  # geom_tippoint(aes(color=source), size=1, shape=20)+
  geom_tiplab(aes(label=tip.lab, color=clade), size=3, offset = .1)+#align = T, linetype = NULL)+
    scale_fill_manual(values=c( "red"), guide='legend')+
  scale_color_manual(values=c("grey30","blue"))+
  theme(legend.position = "bottom")
  
  ggsave(here("figures","ssp_tree.pdf"),p, height=10, width = 12)
  
 
#   # tree to try out node rotations
# p <- 
# ggtree(tree)+
#   geom_tippoint(aes(color=source), size=0.8)+
#   geom_tiplab(aes(label=sp, color=source), size=3, offset = 2)+#align = T, linetype = NULL)+
# geom_nodelab(aes(label=node))+
#   scale_fill_manual(values=c("white", "blue"), guide='legend') 
# p  
# os <- offspring(tree,53)[1:2]
# p <- p%>%flip(os[1],os[2])#%>%rotate(408)%>%rotate(409)%>%rotate(410)
# # p+geom_nodelab(aes(label=node))
# p 

# #####################
# # unrooted trees
p <-
  ggtree(tree, aes(color=clade),size=.5, layout = 'daylight')+
  geom_point2(aes(subset=!isTip ,
                  fill=uf.boot),
              # fill=cut(uf.boot, c(0, 70, 100))),
              shape=21, size=1.5, color="white",stroke = 0) +
  # geom_tippoint(aes(color=source), size=1, shape=20)+
  geom_tiplab(aes(label=tip.lab, color=clade), size=3, offset = .1)+#align = T, linetype = NULL)+
  scale_fill_manual(values=c( "red"), guide='legend')+
  scale_color_manual(values=c("grey","blue"))
# p <- p%>%rotate(405)%>%rotate(408)%>%rotate(409)%>%rotate(410)
ggsave(here("figures","unrooted_ssp_tree.pdf"),p, height=10, width = 7)
# 
# 
# #####################
# # circular trees
# # root.t <- rootnode(tree)
# 
p <-
  ggtree(tree, aes(color=clade),size=.5, layout = 'circular',branch.length = "none")+
  geom_point2(aes(subset=!isTip ,
                  fill=uf.boot),
              # fill=cut(uf.boot, c(0, 70, 100))),
              shape=21, size=1.5, color="white",stroke = 0) +
  # geom_tippoint(aes(color=source), size=1, shape=20)+
  geom_tiplab(aes(label=tip.lab, color=clade), size=3, offset = 1)+#align = T, linetype = NULL)+
  scale_fill_manual(values=c( "red"), guide='legend')+
  scale_color_manual(values=c("grey","blue"))
  ggsave(here("figures","ssp_circular_tree.pdf"),p, height=13, width = 13)
# 
# 
#   ggtree(tree, aes(color=clade),size=.5, layout = 'circular',branch.length = "none")+
#   geom_tippoint(aes(color=source), size=0.1)+
#   geom_tiplab( aes(label=tip.lab,color=source),size=.2)+
#   ggsave(here("circular_tree_labels.pdf"))
