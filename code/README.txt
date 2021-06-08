
(1) get_faa.R
import full genome proteins from NCBI assemblies of phage (infecting Bacillus and Clostridium) and bacteria (based on Burton et al. list) 

(2) hsearch.R
Search all genomic proteins for the major domains of sigma factors region 2 (r2) and region 4 (r4) using HMMER's hmmsearch. The queries used are PFAM HMMs of r2  (PF04542.15) and two models of r4 (PF04545.17, PF08281.13). We then keep only proteins that have hits to both regions (r2 and at least one of the r4)
Filtering based on Paget, M. S. (2015). Bacterial sigma factors and anti-sigma factors: structure, function and distribution. Biomolecules, 5(3), 1245-1265.
"Despite their large variation in size, from ~70 kDa for Group 1 to ~20 kDa for Group 4, all members of the σ70 family possess the σ2 and σ4 domains that include the major RNAP- and promoter-binding determinants "

(3) fasta_headers.R
Check for duplicates. Make minimal and uniform fasta headers for alignmenent and tree construction and save fasta (sigmas_to_align.faa). Save otther data as csv file (sigmas_to_align.csv).


(4) align_trim-tree.sh
align with MAFFT (einsi) and trim alignment with trimAL (automated).
Model selection and ML tree with IQtree2 (default)

# got a warning from iqtree2:
# Number of parameters (K, model parameters and branch lengths): 880
# Sample size (n, alignment length): 141
# Given that K>=n, the parameter estimates might be inaccurate. Thus, phylogenetic estimates should be interpreted with caution.
# based on reponses in IQtree google forum (below) do multiple  (50x) ML tree runs and select best ML (-reruns option does this)
# https://groups.google.com/g/iqtree/c/l8Pi_Xe-Q5A/m/TCNR_mvIAAAJ
# https://groups.google.com/g/iqtree/c/uGeqBo2xm0c/m/BCkAFH46AQAJ
Implemented in code/batch-IQmulti.sh called from align-trim-tree.sh

(5) get_bacterial_features.R
To facilitate tree plotting and analysis I obtained the species specific feature annotations for the aligned proteins. (the sequences downloaded in "get_faa.R" are annotated in MULTISPECIES scheme omplemented by NCBI)

(6) tree-split.R
To further address the "K>=n" warning of iqtree I will implement one of the suggestions in the warning message,  to "Remove the least important sequences from the alignment". There is a large group of ECF sigmas which my previuos experience suggests do not contain phage sigmas (from the set of phages analyzed here). In this R script I identify and filter out the ECF group of sigma factors, as weel as the rpoD (=sigA) house keeping sigmas. From these groups I retained the B. subtilis ECFs and sigA for reference. The filtered set of sequences, along with filtered meta data are written to "data/reduced_set_to_align".

(7) 
rerun alignment, trimming and model selsction as in (4). 
Still grtting "K>n" warning, rhough K is much closer to n now. Nevertheless, I did 50 runs of IQtree to find best ML tree.

(8) 
Best ML tree was the used as basis for UFbootstrap and alRT-SH tests of support for tree.

