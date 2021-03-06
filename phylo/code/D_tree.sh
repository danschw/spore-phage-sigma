#--------------------------
# phylogeny with  RAxML-NG
#--------------------------


#Interactive job on Carbonate
srun -p interactive -N 1 --ntasks-per-node=1 --cpus-per-task=8 --time=07:59:00 --pty bash

#### load dependencies ####

module load raxmlng
  # Flex 2.6.4 loaded.
  # CMake version 3.8.0 loaded.
  # raxmlng version 0.9.0-pthreads loaded.

##### Define paths #####
PARENT=~/GitHub/spore-phage-sigma/phylo
ODIR=${PARENT}/data/align-trim-tree/check_msa
mkdir -p $ODIR
ALN=$PARENT/data/align-trim-tree/sigmas_MafftEinsi.trim

#modeltest-ng selected model : LG+G4
cd $ODIR

raxml-ng --check --msa $ALN  \
--model LG+G4 --data-type AA \
--prefix check-msa 

# WARNING: Duplicate sequences found: 35
# 
# NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
# NOTE: was saved to: /geode2/home/u020/danschw/Carbonate/GitHub/spore-phage-sigma/phylo/data/align-trim-tree/check_msa/check-msa.raxml.reduced.phy

# checked that sequenced removed as duplicates do not include sequences directly discussed in MS
# duplicated-seqs.R

ALN=$ODIR/check-msa.raxml.reduced.phy

# For large alignments, we recommend using the --parse command after, or, instead of
# --check:
raxml-ng --parse --msa $ALN --data-type AA --model LG+G4 --prefix parse-msa 

# ALN=$ODIR//parse-msa.raxml.rba

# * Estimated memory requirements                : 81 MB
# 
# * Recommended number of threads / MPI processes: 3

# get tree
ODIR=$PARENT/data/align-trim-tree/tree
mkdir -p $ODIR
cd $ODIR
raxml-ng --msa $ALN --model LG+G4 --threads 3 --seed 123

mv ../check_msa/check-msa.raxml.reduced.phy.raxml* .
