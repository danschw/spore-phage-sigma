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
ODIR=${PARENT}/data/align-trim-tree/iqtree1
mkdir -p ${ODIR}
cd ${ODIR}
cp ${PARENT}/data/align-trim-tree/sigmas_MafftEinsi.trim ${ODIR}
$IQT -s sigmas_MafftEinsi.trim --seqtype AA -T AUTO

rm sigmas_MafftEinsi.trim 
