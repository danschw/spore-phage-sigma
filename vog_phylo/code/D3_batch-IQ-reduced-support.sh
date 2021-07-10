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
ODIR=${PARENT}/data/reduced_set_to_align/iqtree1-support

mkdir -p ${ODIR}
cd ${ODIR}
cp ${PARENT}/data/reduced_set_to_align/sigmas_MafftEinsi.trim ${ODIR}
cp ${PARENT}/data/reduced_set_to_align/multi-run-iqtree/sigmas_MafftEinsi.trim.treefile ${ODIR}


$IQT -s sigmas_MafftEinsi.trim -m LG+G4  -t sigmas_MafftEinsi.trim.treefile \
--seqtype AA -T AUTO --ufboot 10000 -alrt 10000 --bnni 


#clean up
rm sigmas_MafftEinsi.trim 

