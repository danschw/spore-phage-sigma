#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-26:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=multi-IQ

######  Module commands #####
TOOLS=/N/u/danschw/Carbonate/my_tools
IQT=${TOOLS}/iqtree-2.1.3-Linux/bin/iqtree2

######  Job commands go below this line #####
PARENT=/N/u/danschw/Carbonate/GitHub/spore-phage-sigma/vog_phylo
ODIR=${PARENT}/data/align-trim-tree/multi-run-iqtree1
mkdir ${ODIR}
cd ${ODIR}

cp ${PARENT}/data/align-trim-tree/sigmas_MafftEinsi.trim ${ODIR}

$IQT -s sigmas_MafftEinsi.trim -m Q.pfam+R6 --runs 50 --seqtype AA -T AUTO

rm sigmas_MafftEinsi.trim 
