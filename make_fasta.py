#!/usr/bin/python

# Turn table of effector sequences (Table S3 Thilliez et al. 2018)
# Into fasta format

fh = open("capsici_effectors.txt", "r")

for line in fh:
	line = line.strip().split("\t")
	name = line[0]
	seq = line[3]
	print ">" + name
	print seq
fh.close()
