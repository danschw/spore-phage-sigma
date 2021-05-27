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
  for (f in bacteria.faa[20]){
    
    faa.name <- 
      str_remove(f,".faa.gz")%>%
      str_remove(".*/")
    
    hmm.name <- 
      str_remove(h,".hmm$")%>%
      str_remove(".*/")
    
    filename <- paste0("hsearch_results/bacteria/",hmm.name,"_X_",faa.name,".txt")
    
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

hits2 <- hits %>%
  select(target.name, query.name) %>%
  mutate(present = 1) %>% # create a dummy column
  pivot_wider(names_from = query.name, values_from = present, values_fill = 0 ) %>%
  mutate(both = Sigma70_r2==1 & (Sigma70_r4==1|Sigma70_r4_2==1))
#--------------
#--------------
# need to filter hits: containing both r2 and r4 domains
#--------------
#--------------


write_csv(hits,here("hmm_hits.csv"))

# download AA fasta files of hits from NCBI protein data base
for (i in unique(hits$target.name)){
  #get the fasta using linux
  wsl <- paste0("wsl efetch -db protein -id ",i," -format fasta")
  fa <- shell(wsl,intern=TRUE)
  #wirte to a file
  write(fa, file = here("hmm_hits.faa"), append = TRUE)
  print(i)
  Sys.sleep(1) #to avoid limits by NCBI
}

# 429 Too Many Requests
# PLEASE REQUEST AN API_KEY FROM NCBI
# No do_post output returned from 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=AMO25912.1&rettype=fasta&retmode=text&edirect_os=linux&edirect=12.0&tool=edirect&email=danschw@bl-bio-g33qhl2.localdomain'
# Result of do_post http request is
