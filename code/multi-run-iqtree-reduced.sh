#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=0-23:59:00
#SBATCH --mem=50gb
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=multi-IQ

######  Module commands #####
TOOLS=/N/u/danschw/Carbonate/my_tools
IQT=${TOOLS}/iqtree-2.1.3-Linux/bin/iqtree2

######  Job commands go below this line #####
PARENT=/N/u/danschw/Carbonate/GitHub/spore-phage-sigma
ODIR=${PARENT}/data/reduced_set_to_align/multi-run-iqtree
mkdir ${ODIR}
cd ${ODIR}

cp ${PARENT}/data/reduced_set_to_align/sigmas_MafftEinsi.trim ${ODIR}

$IQT -s sigmas_MafftEinsi.trim -m LG+G4 --runs 50 --seqtype AA -T AUTO

#$IQT -s sigmas_MafftEinsi.trim --runs 50 --seqtype AA -T AUTO

rm sigmas_MafftEinsi.trim 

