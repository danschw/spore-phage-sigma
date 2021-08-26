#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:30:00
#SBATCH --mem=1gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=rax-bs

##### load dependencies #####
module load raxmlng

##### Define paths #####
PARENT=~/GitHub/spore-phage-sigma/phylo
ODIR=${PARENT}/data/reduced_set_to_align/tree/bootstraps
mkdir -p $ODIR
ALN=$PARENT/data/reduced_set_to_align/check_msa/parse.raxml.rba


##### run raxml-ng #####
SEED=$1$1$1

cd $ODIR

raxml-ng --bootstrap --msa $ALN --model LG+G4  --prefix bs$SEED --seed $SEED --threads 2 \
--tree pars{5} --bs-trees 25