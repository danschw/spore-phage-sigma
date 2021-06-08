#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=0-6:59:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=IQ-thread

######  Module commands #####
TOOLS=/N/u/danschw/Carbonate/my_tools
IQT=${TOOLS}/iqtree-2.1.3-Linux/bin/iqtree2


######  Job commands go below this line #####
PARENT=/N/u/danschw/Carbonate/GitHub/spore-phage-sigma
ODIR=${PARENT}/data/reduced_set_to_align/iqtree-support

mkdir -p ${ODIR}
cd ${ODIR}
cp ${PARENT}/data/reduced_set_to_align/sigmas_MafftEinsi.trim ${ODIR}
cp ${PARENT}/data/reduced_set_to_align/multi-run-iqtree/sigmas_MafftEinsi.trim.treefile ${ODIR}


$IQT -s sigmas_MafftEinsi.trim -m LG+G4 -t sigmas_MafftEinsi.trim.treefile --seqtype AA -T AUTO \
--ufboot 10000 -alrt 10000 --bnni

#clean up
rm sigmas_MafftEinsi.trim 
rm sigmas_MafftEinsi.trim.treefile
