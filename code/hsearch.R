library(here)
library(tidyverse)

# this script relies on the windows subsytem for Linux (wsl)
# on which hmmer has been installed
# paths passed to wsl system cannot be constructed by here()
setwd(here("data"))


# list of sigma factor domain hmmms from PFAM
hmm <- dir("sigma_hmms/", full.names = T)

#-------------------------#
# search in phage genomes #
#-------------------------#
#proteins to be searched 
phage.faa <- dir("phage_faa/" ,pattern = ".faa.gz$",full.names = T)

# make directory for results
if (!dir.exists(here("data","hsearch_results/phage"))){
  dir.create(here("data","hsearch_results/phage"), recursive = T)
} 


for (h in hmm){
  for (f in phage.faa){
    
    faa.name <- 
      str_remove(f,".faa.gz")%>%
      str_remove(".*/")
    
    hmm.name <- 
      str_remove(h,".hmm$")%>%
      str_remove(".*/")
    
    filename <- paste0("hsearch_results/phage/",hmm.name,"_X_",faa.name,".txt")
    
    #PFAM hmms have non-numeric name and contain noise cutoff (--cut_nc)
    if (is.na(as.numeric(hmm.name))){
      shell(paste("wsl hmmsearch --noali --cut_nc --tblout", filename, h, f))
      next
    } 
    # superfam hmms have numeric names and no NC noise cutoff
    shell(paste("wsl hmmsearch --noali -T 20 --tblout", filename, h, f))

  }
}


#-----------------------------#
# search in bacterial genomes #
#-----------------------------#
#proteins to be searched 
bacteria.faa <- dir("bacteria_faa/" ,pattern = ".faa.gz$",full.names = T)

# make directory for results
if (!dir.exists(here("data","hsearch_results/bacteria"))){
  dir.create(here("data","hsearch_results/bacteria"), recursive = T)
} 


for (h in hmm){
  for (f in bacteria.faa){
    
    faa.name <- 
      str_remove(f,".faa.gz")%>%
      str_remove(".*/")
    
    hmm.name <- 
      str_remove(h,".hmm$")%>%
      str_remove(".*/")
    
    filename <- paste0("hsearch_results/bacteria/",hmm.name,"_X_",faa.name,".txt")
    
    #PFAM hmms have non-numeric name and contain noise cutoff (--cut_nc)
    if (is.na(as.numeric(hmm.name))){
      shell(paste("wsl hmmsearch --noali --cut_nc --tblout", filename, h, f))
      next
    } 
      # superfam hmms have numeric names and no NC noise cutoff
      shell(paste("wsl hmmsearch --noali -T 20 --tblout", filename, h, f))
    
  }
}

#--------------------
# phage
#------------------
# parse hmmserch results
l.res <- list.files(here("data/hsearch_results/phage"), full.names = T)
# variable to store reults
agregate <- character()

for (i in 1:length(l.res)){
  current <- readLines(l.res[i])#
  empty.line <- which(current=="#")
  # If no results move to next file
  if(empty.line ==4) next
  
  tmp <- read_lines(l.res[i],skip = 3,n_max = empty.line-4)%>%
    str_replace_all("  * ",",")%>%
    # there is one space that is missed
    str_replace(" ",",")
  
  agregate <- c(agregate,tmp)
}

write(agregate,here("data","phage_hmm_hits.csv"))

#--------------------
# bacteria
#------------------
# parse hmmserch results
l.res <- list.files(here("data/hsearch_results/bacteria"), full.names = T)
# variable to store reults
agregate <- character()

for (i in 1:length(l.res)){
  current <- readLines(l.res[i])#
  empty.line <- which(current=="#")
  # If no results move to next file
  if(empty.line ==4) next
  
  tmp <- read_lines(l.res[i],skip = 3,n_max = empty.line-4)%>%
    str_replace_all("  * ",",")%>%
    # there is one space that is missed
    str_replace(" ",",")
  
  agregate <- c(agregate,tmp)
}

write(agregate,here("data","bacteria_hmm_hits.csv"))
#--------------

cnames <- c("target.name","target.accession","query.name","query.accession","full.seq.E-value","full.seq.score","full.seq.bias","best.domain.E-value","best.domain.score","best.domain.bias","exp","reg","clu","ov","env","dom","rep","inc","description.of.target")

hits.p <- read_csv(here("data","phage_hmm_hits.csv"), col_names = F)
hits.b <- read_csv(here("data","bacteria_hmm_hits.csv"), col_names = F)

colnames(hits.p) <- cnames
colnames(hits.b) <- cnames

hits.b$group <- "bacteria"
hits.p$group <- "phage"

hits <- bind_rows(hits.b, hits.p)

# save hits table
write_csv(hits,here("data/hmm_ALL_hits.csv"))
#hits <- read_csv(here("data/hmm_ALL_hits.csv"))



#--------------
# need to filter hits: containing both r2 and r4 domains
# hits2 <- hits %>%
#   select(target.name, query.name, group) %>%
#   mutate(present = 1) %>% # create a dummy column
#   pivot_wider(names_from = query.name, values_from = present, values_fill = 0 ) %>%
#   mutate(both = Sigma70_r2==1 & (Sigma70_r4==1|Sigma70_r4_2==1))

# names of superfamily HMMs used
super.name <- read_tsv(here("data/superfam_sigma_names.tsv"),
                       col_names = F)
super.name <- super.name %>% 
  select(query.name = 1, name = 5) %>% 
  mutate(domain = case_when(str_detect(name, "2") ~ "r2_sf",
                            str_detect(name, "4") ~ "r4_sf")) %>% 
  select(-name)

hits2 <- hits %>%
    select(target.name, description.of.target, query.name, group) %>% 
    left_join(., super.name ) %>% 
    mutate(domain = case_when(!is.na(domain) ~ domain,
                              str_detect(query.name, "r2") ~ "r2_pfam",
                              str_detect(query.name, "r4") ~ "r4_pfam")) %>% 
  mutate(present = 1) %>% # create a dummy column
  select(-query.name, -group) %>% 
  distinct() %>% 
  pivot_wider(names_from = domain, values_from = present, values_fill = 0 ) %>% 
  mutate(r2 = ((r2_pfam + r2_sf) > 0), r4 = ((r4_pfam + r4_sf) > 0)) %>% 
  mutate(both = r2 & r4)





keep.hits <- hits2 %>% 
  filter(both==TRUE) %>% 
  pull(target.name)

hits.filt <- hits %>% 
  filter(target.name %in% keep.hits)

# save filtered hits table
write_csv(hits.filt,here("data/hmm_r2-r4_hits.csv"))
#--------------





# download AA fasta files of hits from NCBI protein data base
if (file.exists(here("data/hmm_sigmas_hits.faa"))){
  file.remove(here("data/hmm_sigmas_hits.faa"))
}

for (i in unique(hits.filt$target.name)){
  #get the fasta using linux
  wsl <- paste0("wsl efetch -db protein -id ",i," -format fasta")
  fa <- shell(wsl,intern=TRUE)
  #wirte to a file
  write(fa, file = here("data/hmm_sigmas_hits.faa"), append = TRUE)
  print(i)
  Sys.sleep(2) #to avoid limits by NCBI
}

#cleanup
file.remove(here("data/bacteria_hmm_hits.csv"))
file.remove(here("data/phage_hmm_hits.csv"))
file.remove(here("data/tmp_hmm_hits.csv"))



# ##########################
# hits2 <- hits2 %>% 
#   mutate(old.keep = (r2_pfam + r4_pfam ) == 2) %>% 
#   mutate(sig_string = str_detect(description.of.target, regex("sigma", ignore_case = T)))
