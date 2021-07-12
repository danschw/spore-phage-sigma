#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-6:59:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=IQ-support

######  Module commands #####
TOOLS=/N/u/danschw/Carbonate/my_tools
IQT=${TOOLS}/iqtree-2.1.3-Linux/bin/iqtree2


######  Job commands go below this line #####
PARENT=/N/u/danschw/Carbonate/GitHub/spore-phage-sigma/vog_phylo
ODIR=${PARENT}/data/curated_set_to_align/iqtree-fullALN-support

mkdir -p ${ODIR}
cd ${ODIR}
cp ${PARENT}/data/curated_set_to_align/sigmas_MafftEinsi.aln ${ODIR}


$IQT -s sigmas_MafftEinsi.aln -m Q.pfam+F+R4  \
--seqtype AA -T AUTO --ufboot 10000 -alrt 10000 --bnni 


#clean up
rm sigmas_MafftEinsi.aln 

