#!/bin/bash

#This was executed on Carbonate interactive job


#### local dependencies ####
TOOLS=/N/u/danschw/Carbonate/my_tools
MAFFT=${TOOLS}/mafft-7.475-without-extensions/bin
# MAFFT v7.475 (2020/Nov/23)
# MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
# https://mafft.cbrc.jp/alignment/software/

TRIMAL=/geode2/home/u020/danschw/Carbonate/my_tools/trimal-trimAl/source
# trimAl v1.4.rev22 build[2015-05-21]. 2009-2015. Salvador Capella-Gutierrez and Toni Gabaldón.
# trimAl webpage: http://trimal.cgenomics.org
# (used in Consuelo Gazitúa wt (2020) paper (Sullivan lab))


##### Define paths #####
PARENT=~/GitHub/spore-phage-sigma
ODIR=${PARENT}/data/align-trim-tree
mkdir -p ${ODIR}


# alignment
cd $MAFFT

./mafft-einsi  ${PARENT}/data/sigmas_to_align.faa > ${ODIR}/sigmas_MafftEinsi.aln

# trim alignment
cd $TRIMAL

./trimal -in ${ODIR}/sigmas_MafftEinsi.aln -out ${ODIR}/sigmas_MafftEinsi.trim -automated1



# local iqtree2 with model search
#IQT=/N/u/danschw/Carbonate/my_tools/iqtree-2.1.3-Linux/bin/iqtree2
# # IQ-TREE multicore version 2.1.3 COVID-edition for Linux 64-bit built Apr 21 2021
#cd $ODIR
#$IQT -s sigmas_MafftEinsi.trim -T auto

sbatch ${PARENT}/code/batch-IQ.sh 
# Best-fit model: Q.pfam+R6 chosen according to BIC


# got a warning in sigmas_MafftEinsi.trim.iqtree:
# Number of parameters (K, model parameters and branch lengths): 939
# Sample size (n, alignment length): 91
# Given that K>=n, the parameter estimates might be inaccurate. Thus, phylogenetic estimates should be interpreted with caution.
# based on reponses in IQtree google forum (links below) doing multipe  ML tree runs.
# https://groups.google.com/g/iqtree/c/l8Pi_Xe-Q5A/m/TCNR_mvIAAAJ
# https://groups.google.com/g/iqtree/c/uGeqBo2xm0c/m/BCkAFH46AQAJ

sbatch ${PARENT}/code/multi-run-iqtree.sh



