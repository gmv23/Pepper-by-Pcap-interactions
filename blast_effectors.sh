#!/usr/bin/bash

#Blast list of effector sequences to genome to pull out coordinates of effectors

path_to_genome="/workdir/gmv23/sexual_pop/genome"

blastn -db $path_to_genome/Phyca11_unmasked_genomic_scaffolds \
-query Pc_effectors.fa -out Pc_effectors_blastout.txt -outfmt 6 \
-max_target_seqs 1 -max_hsps 1

