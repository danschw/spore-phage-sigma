#!/bin/bash

#This was executed on Carbonate


#### load dependencies ####
module load muscle
	#muscle version 3.8.31 loaded.


##### Define paths #####
PARENT=~/GitHub/spore-phage-sigma/

muscle -in $PARENT/data/sigmas_to_align.faa -out $PARENT/data/sigmas_autoMuscle.aln
