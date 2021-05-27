library(here)
library(tidyverse)


# genbank viral index file
d.gbv <- read_tsv(here("data","assembly_summary_20201104.txt"),skip = 1)


filt <- 
d.gbv%>%
  filter(str_detect(organism_name,"^Bacillus")|
           str_detect(organism_name,"^Clostridium"))

if (!dir.exists(here("protein_faa/")))
    dir.create(here("protein_faa/"))
setwd(here("protein_faa/"))

for (i in 259:nrow(filt)){
  # problem with phage GCA_002609465.1 (i=90) Bacillus phage Crookii	N
  if (filt$assembly_accession[i]=="GCA_002609465.1") next
  # problem with phage GCA_002609465.1 (i=247) Bacillus phage SBSphiJ
  if (filt$assembly_accession[i]=="GCA_002990895.1") next
  # problem with phage GCA_002609465.1 (i=248) Bacillus phage SBSphiC
  if (filt$assembly_accession[i]=="GCA_002990925.1") next
  # problem with phage GCA_002609465.1 (i=258) Bacillus phage EZ-2018a	
  if (filt$assembly_accession[i]=="GCA_003012515.1") next
  file.name <- paste0(str_remove(filt$ftp_path[i],".*/"),"_protein.faa.gz")
  
  ftp <- paste0(filt$ftp_path[i],"/",file.name)
  
  download.file(ftp, destfile = paste0(filt$organism_name[i],".faa.gz"))
}


#fix names to remove spaces
for (n in list.files()){
  file.rename(n,str_replace_all(n," ","_"))
}

####################################################
# Get bacterial faa files specified by Burton et al. 2019
d.burton <- read_csv(here("Burton_S6.csv"))

#fix names to remove spaces
d.burton <- d.burton%>%
  mutate(organism_name=str_replace_all(organism_name," ","_"))

setwd(here("bacteria_faa/"))

for (i in 2:nrow(d.burton)){

  file.name <- paste0(str_remove(d.burton$ftp_path[i],".*/"),"_protein.faa.gz")
  
  ftp <- paste0(d.burton$ftp_path[i],"/",file.name)
  
  download.file(ftp, destfile = paste0(d.burton$organism_name[i],".faa.gz"))
}

