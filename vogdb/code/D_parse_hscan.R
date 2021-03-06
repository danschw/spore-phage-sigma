library(here)
library(tidyverse)
library(cowplot)
library(scales)
source(here("vogdb/code/parse_hmmer_tbl.R"))

# parse hmmserch results
l.res <- list.files(here("vogdb/data/hscan_vogXtigr/"), full.names = T)
# variable to store reults
hits.all <- tibble()


for (i in 1:length(l.res)){
  cur.hits <- read_tblout(l.res[i])
  
  # make row for 0 hits
  if(nrow(cur.hits)==0){
    #add row
    cur.hits[1,] <- NA
    # get protein name
    cur.hits$query_name[1] <- 
      str_remove(l.res[i],".*hscan_vogXtigr/")%>%
      str_remove(".txt")
    #assign no_hit fam
    cur.hits$description[1] <- "no_hit:no_hit"
    
    #assign non-significant values to E-values
    cur.hits <- cur.hits%>%
      mutate(across(where(is.numeric),~1))
  }
  
  hits.all <- bind_rows(hits.all,cur.hits)
  

}

# separate TIGR descriptions
hits.all <- hits.all%>%
  separate(description,into = c("id","description"),sep = ":" )


# need to solve the no hit rows which do not have taxon
# #separate query name
#   separate(query_name, 
#            into = c("taxon","protein"),
#            sep = "\\.",
#            extra="merge", remove = FALSE)


##############
# get 3 best hits for each protein
# defining best by sequence evalue()

best3 <- hits.all%>%
  group_by(query_name)%>%
  mutate(sequence_rate=-log10(sequence_evalue))%>%
  slice_max(sequence_rate, n=3)%>%
  # define hit order
  arrange(desc(sequence_rate))%>%
  mutate(place=row_number())%>%
  mutate(place=paste0("tigr.hit",place))%>%
  ungroup()

# hist(best3$sequence_rate, breaks = 100)
# all(!is.na(best3$sequence_rate))
# all(!is.na(best3$sequence_rate))

# find the top most specific hit
# if best hit is a genral term:
# "sigma70-ECF"
                # "SigBFG"
# look down the list to find a more specific one.
# if a more specific one does is not present use the general one
# in case the 2 general terms are listed and no 3rd term than use "SigBFG" as it is more specific


best3<- best3%>%
  pivot_wider(query_name, names_from=place,values_from=id)%>%
  add_column(tigr.hit=NA)%>%
  mutate(tigr.hit=as.character(tigr.hit))



g.terms <- c("sigma70-ECF","SigBFG")
# g.terms <- c("sigma70-ECF")

for(i in 1:nrow(best3)){
  if(!best3$tigr.hit1[i] %in% g.terms|
     is.na(best3$tigr.hit2[i])){
    best3$tigr.hit[i] <- best3$tigr.hit1[i]
    next
  }
  if(!best3$tigr.hit2[i] %in% g.terms|
     is.na(best3$tigr.hit3[i])){
    best3$tigr.hit[i] <- best3$tigr.hit2[i]
    next
  }
  best3$tigr.hit[i] <- best3$tigr.hit3[i]

}

# these rules did the job. 
# SigBFG was completely excluded (always had a more specific term to take over)
# sigma70-ECF was retained only for proteins that had no other hit

best3%>%
  group_by(query_name)%>%
  group_by(tigr.hit)%>%
  summarise(n=n())%>%
  mutate(tigr.hit=fct_reorder(tigr.hit,n))%>%
  ggplot(aes(tigr.hit,n))+
  geom_col()+
  geom_label(aes(label=paste0("n=",n)), y=150)+
  coord_flip()+
  theme_cowplot()+
  xlab("best specific TIGRFAM match")  
  ggsave2(filename = here("vogdb/figures/vogXtigr.png"),
          width = 7,height = 7)


#looking at e-values of best hits
best3.evalue <- hits.all%>%
  group_by(query_name)%>%
  slice_max(sequence_evalue, n=3)%>%
  # define hit order
  arrange(desc(sequence_evalue))%>%
  mutate(place=row_number())%>%
  ungroup()


best3.evalue%>%
  filter(id!="no_hit")%>%
  # pivot_wider(query_name, names_from=place,values_from=sequence_evalue)%>%
  ggplot(aes(x=place, y=-log10(sequence_evalue)))+
  geom_jitter(aes(color=sequence_score),alpha=0.7, width = 0.3)+
  geom_hline(yintercept = c(3), color="red")+
  geom_hline(yintercept = c(5), color="pink")+
  facet_wrap(~id)+
  scale_y_log10()+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")+
  theme_cowplot()


# all have an e-value lower than 1e-3 expect for a single sigH hit
# I don't think there is a problem




################################################
# combining vogdb and virus-host db info on proteins

# load in the data made in previous script
load(here("vogdb","data","vog_sigma_clean_Whost.RData"))

#combining of data frames will be by protein number.
# extract it from hscan query_name
# remove leading taxon number from query name if present
#leaving protein id that starts with letter

best3 <- best3%>%
  mutate(protein=str_remove(best3$query_name, regex("^[0-9]*\\.?")))

##QC
# anyDuplicated(best3$protein) # 0
# is.na(best3$protein)%>%sum() # 0
# all(d.faa$protein %in% best3$protein) #TRUE
# all(best3$protein %in% d.faa$protein) #TRUE


d.faa <- 
  best3%>%
  select(-query_name)%>%
  full_join(d.faa, ., by="protein")


############################################
# Do tigr hits correspond to VOG?
############################################
d.faa%>%
  select(vog,tigr.hit,tigr.hit1,tigr.hit2,tigr.hit3 )%>%
  pivot_longer(-vog, names_to="hit", values_to="tigr.fam")%>%
  group_by(vog,hit,tigr.fam)%>%
  summarise(n=n())%>%
  filter(!is.na(tigr.fam))%>%
  ggplot(aes(x=vog, y=tigr.fam))+
  geom_tile(aes(fill=n))+
  facet_wrap(~hit)+
  scale_fill_viridis_b()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=30,hjust = 1),
        legend.position = "bottom")+
  ggsave2(filename = here("vogdb/figures/vogXtigr_HitPlace.png"),
          width = 10,height =10)

# mostly no clear 1 to 1 correspondence

############################################
# parsing sigma factor type by 
# number of sigma factors in phage
###########################################

# add column of sigma factor content per genome
n.sig <- d.faa%>%
  group_by(taxon, `virus name`)%>%
  summarise(n.sig=n())%>%
  ungroup()

d.faa <- n.sig%>%
  select(-`virus name`)%>%
  right_join(.,d.faa, by="taxon")

###################################
# classes of sigma factors
sig.class <- read_csv(here("vogdb/data/tigr_Bs_sigma.csv"))

d.faa <- sig.class%>%
  select(ID,sig.class, sig.group)%>%
  left_join(d.faa,., by=c("tigr.hit"="ID"))

d.faa$sig.class[d.faa$tigr.hit=="no_hit"] <- "no_hit"
d.faa$sig.group[d.faa$tigr.hit=="no_hit"] <- "no_hit"




##########################
# Plot by phylum
##########################
n.labs <- d.faa%>%
  group_by(protein,phylum)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(phylum)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))

tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")


p.phylum <- 
  d.faa%>%
  # arrange by frequency of sigma factors
  mutate(tigr.hit=fct_infreq(tigr.hit))%>%
  group_by(phylum,tigr.hit)%>%
  summarise(n=n())%>%
    group_by(phylum)%>%
    mutate(perc=n/sum(n))%>%
    ungroup()%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  #add color to sporulation sigma factors
  mutate(spore.sigma=if_else(tigr.hit%in%tigr.spore,"spore","other"))%>%
    ggplot( aes(fct_rev(tigr.hit),perc, group = phylum)) + 
  #line for 0
  geom_hline(yintercept = 0, color="grey")+
  geom_col(aes(fill=spore.sigma)) + 
  geom_text(data=n.labs, aes(label=lab), x=2,y=0.5)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage sigma factor genes") +
  xlab("Sigma Factor TIGRFAM")+
  facet_grid(~phylum)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  coord_flip()+
  panel_border()

p <- p.phylum +
  ggtitle("phage sigma factor type by host phylum")
  ggsave2(here("vogdb","figures","tigr_sigma_HostPhylum.png"),
          plot = p, width = 12,height = 8)

  ###
  n.gene.phylum <- d.faa%>%
    mutate(phylum=str_remove(phylum,"/.*"))%>%
    group_by(phylum)%>%
    summarise(n=n(),.groups = "drop") %>% 
    mutate(phyl.lab = paste0(phylum,"\n(n=",n,")"))
# summarize by spore function
  d.plot <- 
    d.faa%>%
    #classify by spore function
    mutate(sigma = if_else(tigr.hit %in% tigr.spore,tigr.hit, "other" )) %>% 
    mutate(sigma = str_remove(sigma, "spore_")) %>% 
    mutate(sigma = str_replace(sigma, "sigma", "sig")) %>% 
    mutate(sigma = str_replace(sigma, "Sig", "sig")) %>% 
    mutate(sigma = fct_relevel(sigma, "sigF", "sigG", "sigK", "sigE", "other") %>% 
             fct_rev()) %>%
    group_by(phylum,sigma)%>%
    summarise(n=n(), .groups = "drop") %>% 
    group_by(phylum)%>%
    mutate(perc=n/sum(n)) %>% 
    
    mutate(phylum=str_remove(phylum,"/.*"))%>%
    mutate(phylum=fct_reorder(phylum,n)) %>% 
    left_join(., n.gene.phylum, by = "phylum")
  
  p.phylum <- d.plot %>% 
    ggplot(aes(x=phyl.lab, y = perc,fill = sigma)) + 
    geom_bar(position="fill", stat="identity", color = "black", width = 0.5) +
    
    xlab("Host Phylum") +
    ylab(NULL)+
    guides(fill = guide_legend("Sigma\nFactor\nType", reverse = T))+
    
    scale_y_continuous(labels=scales::percent, position = "right") +
    scale_fill_viridis_d(direction = -1) + 
    theme_classic(base_size = 13)+
    panel_border(color = "black")+
    coord_flip()
  
  
    ggsave2(here("vogdb","figures","tigr_Nsigma_HostPhylum.png"),
            p.phylum, width = 4,height = 3)   
    
  



###################################
# Focus on phages with 3 sigma factors
d.nsig <- 
d.faa%>%
  #classify by spore function
  mutate(sigma = if_else(tigr.hit %in% tigr.spore,tigr.hit, "other" )) %>% 
  mutate(sigma = str_remove(sigma, "spore_")) %>% 
  mutate(sigma = str_replace(sigma, "sigma", "sig")) %>% 
  mutate(sigma = str_replace(sigma, "Sig", "sig")) %>% 
  mutate(sigma = fct_relevel(sigma, "sigF", "sigG", "sigK", "sigE", "other") %>% 
           fct_rev()) %>% 
  # assign a positional index to each sigma factor
  group_by(taxon)%>%
  arrange(tigr.hit)%>%
  mutate(sig.position=row_number())%>%
  ungroup()%>%
  mutate(`virus name`=str_remove(`virus name`, ".*phage "))%>%
  mutate(`virus name`=str_remove(`virus name`, ".*virus"))%>%
  mutate(`virus name`=as_factor(`virus name`)) %>% 
      
  separate(family.etc, into = c("family","genus","species","strain"),sep=";") %>% 
  # genus is also in viral sp name as first word
  mutate(genus2=str_extract(sp,regex(".*? ")) %>% trimws()) %>% 
  mutate(genus = trimws(genus)) %>% 
  mutate(genus.plot = if_else(is.na(genus), genus2, genus))
  
    
  #plot
  pA <-   d.nsig %>% 
      filter(genus.plot == "Bacillus") %>% 
  # arrange phage by n.sigmas
  mutate(`virus name` = fct_reorder(`virus name`, n.sig)) %>% 
  ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
  geom_tile(aes(fill=sigma), color="black")+
  theme_classic()+
  panel_border(color = "black")+
  scale_fill_viridis_d(direction = -1)+
  facet_grid(genus.plot~., scales = "free", space = "free")+
    geom_text(aes(label = genus.plot), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title.y = element_blank())+
    xlab("Sigma Factor Gene")+
    guides(fill = guide_legend("Sigma\nFactor\nType", reverse = T))
  

  
  pB <-   d.nsig %>% 
    mutate(genus.plot = fct_infreq(genus.plot)) %>% 
    filter(phylum == "Firmicutes") %>% 
    filter(genus.plot != "Bacillus") %>% 
   
    # arrange phage by n.sigmas
    mutate(`virus name` = fct_reorder(`virus name`, n.sig)) %>% 
    ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
    geom_tile(aes(fill=sigma), color="black")+
    theme_classic()+
    panel_border(color = "black")+
    scale_fill_viridis_d(direction = -1)+
    facet_grid(genus.plot~., scales = "free", space = "free")+
    geom_text(aes(label = genus.plot), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          legend.position = "none",
          axis.title.y = element_blank())+
    xlab("Sigma Factor Gene")+
    xlim(NA,3)
  

p <- plot_grid(pB, pA, rel_widths = c(1,1.4))  
  
ggsave2(here("vogdb/figures/","sigma_TIGR_content_Firmi.png"),
        plot = p, width = 8, height = 11.5)


# 
# pC <-   d.nsig %>% 
#   filter(phylum != "Firmicutes") %>% 
#   # arrange phage by n.sigmas
#   mutate(`virus name` = fct_reorder(`virus name`, n.sig)) %>% 
#   ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
#   geom_tile(aes(fill=sigma), color="black")+
#   theme_classic()+
#   panel_border(color = "black")+
#   scale_fill_viridis_d(direction = -1)+
#   facet_grid(phylum~., scales = "free", space = "free")+
#   geom_text(aes(label = phylum), x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
#   theme(strip.background = element_blank(),
#         strip.text = element_blank(),
#         axis.title.y = element_blank())+
#   xlab("Sigma Factor Gene")+
#   guides(fill = guide_legend("Sigma\nFactor\nType", reverse = T))
# 
# pC
# Save data ---------------------------------------------------------------
  
  save(d.faa, file = here("vogdb/data/vog_sigma_clean_tigr.RData"))

  write_csv(d.faa %>% select(-seq), file = here("vogdb/data/vog_sigma_clean_tigr.csv"))
# ###############
# d.nsig3%>%
#   # assign a positional index to each sigma factor
#   group_by(taxon)%>%
#   # arrange(sigma.group)%>%
#   mutate(sig.position=row_number())%>%
#   ungroup()%>%
#   mutate(`virus name`=str_replace(`virus name`, "phage", "\u03c6"))%>%
#   mutate(`virus name`=str_replace(`virus name`, "virus", "\u03c6"))%>%
#   mutate(`virus name`=as_factor(`virus name`))%>%
#   #plot
#   ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
#   geom_tile(aes(fill=sigma.group), color="black")+
#   theme_cowplot()+
#   scale_fill_viridis_d()+
#   facet_wrap(~n.sig, scales = "free")+
#   theme(axis.text.x = element_blank())+
#   ggsave2(here("vogdb/figures/","sigma_pfam_content_firmicutes.png"),
#           width = 13, height = 18)




