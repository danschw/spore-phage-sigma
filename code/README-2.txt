
(1) get_faa.R
import full genome proteins from NCBI assemblies of phage (infecting Bacillus and Clostridium) and bacteria (based on Burton et al. list) 

(2) hsearch.R
Search all genomic proteins for the major domains of sigma factors region 2 (r2) and region 4 (r4) using HMMER's hmmsearch. The queries used are PFAM HMMs of r2  (PF04542.15) and two models of r4 (PF04545.17, PF08281.13). We then keep only proteins that have hits to both regions (r2 and at least one of the r4)
Filtering based on Paget, M. S. (2015). Bacterial sigma factors and anti-sigma factors: structure, function and distribution. Biomolecules, 5(3), 1245-1265.
"Despite their large variation in size, from ~70 kDa for Group 1 to ~20 kDa for Group 4, all members of the σ70 family possess the σ2 and σ4 domains that include the major RNAP- and promoter-binding determinants "

(3) fasta_headers.R
Check for duplicates. Make minimal and uniform fasta headers for alignmenent and tree construction and save fasta (sigmas_to_align.faa). Save otther data as csv file (sigmas_to_align.csv).

###############
(4) align-muscle.sh
Align with Muscle (default options) on carbonate (sigmas_autoMuscle.aln)

(5) phylo-iqtree2.sh
###############

(4) align_trim-tree.sh
align with MAFFT (einsi) and trim alignment with trimAL (automated).
model selection and ML tree with IQtree2 (default)

