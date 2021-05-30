library(tidyverse)
library(here)
library(cowplot)
library(seqinr)

# import data of sigmas generated in hmmer search
hits <- read_csv(here("data/hmm_r2-r4_hits.csv"))

# import fasta file of hit protein sequences
sigma_fa <- read.fasta(here("data/hmm_sigmas_hits.faa"),whole.header = TRUE)

# parse fasta file into table
d.faa <- tibble()

for(i in 1:length(sigma_fa)){
  faa <- sigma_fa[i]
  
  d.faa <- 
    tibble(header=getName(faa),seq=getSequence(faa))%>%
    mutate(protein=str_extract(header,".*\\.. "))%>%
    mutate(protein=str_remove(protein, " $"))%>%
    mutate(sp=str_extract(header,"\\[.*\\]"))%>%
    mutate(sp=str_remove(sp, "\\["))%>%
    mutate(sp=str_remove(sp, "\\]"))%>%
    mutate(description=str_extract(header," .*\\["))%>%
    mutate(description=str_remove(description, " \\["))%>%
    bind_rows(d.faa)
}

#add column of bacteria or phage
d.faa <- hits %>% 
  select(target.name, group) %>% 
  distinct() %>% 
  left_join(d.faa, ., by = c("protein"="target.name"))
###############
# removing duplicates
duplicated(d.faa$protein)%>%sum() # 0, no duplicate protein entries
d.faa %>% 
  filter (duplicated(seq)) %>% 
            group_by(group) %>% 
            summarise(n=n())
# 31 duplicate sequences, all in phages

d.faa%>%
  group_by(seq)%>%
  summarise(n=n())%>%
  filter(n>1) %>% 
    arrange(desc(n))
# there are 22 sequences with 2-4 duplicates

d.faa%>%
  group_by(seq, sp)%>%
  summarise(n=n())%>%
  filter(n>1) %>% 
  arrange(desc(n))
# but none of th duplicated sequences are from the same phage, keeping all


#make fasta headers
d.faa <- d.faa %>% 
  mutate(new.header = paste(protein, group, sep = "_"))

write.fasta(file.out = here("data","sigmas_to_align.faa"),
            sequences = d.faa$seq,
            names = d.faa$new.header)

#save data on sequences used for alignment
d.faa %>% 
  select(-seq) %>% 
  write_csv(here("data","sigmas_to_align.csv"))
