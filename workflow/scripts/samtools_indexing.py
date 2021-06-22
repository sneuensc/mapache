#!/usr/bin/env python

import os

## get the reference
fasta=snakemake.input.orig
filename, file_extension = os.path.splitext(fasta)
orig_dir = os.path.abspath(os.path.dirname(filename))
orig_prefix = os.path.basename(filename)

## get the new reference folder
new_prefix = snakemake.wildcards.id_genome
new_dir=os.path.abspath('results/00_reference/' + new_prefix)

## check if the indexes are present
ext = [".fai"] 
files1 = [filename + s for s in ext]
files2 = [fasta + s for s in ext]
if all([os.path.isfile(f) for f in files1]):    ## foo.ext
	for item in files1:
		path_to = os.path.join(new_dir, os.path.basename(item).replace(orig_prefix, new_prefix + '.fasta'))
		os.symlink(item, path_to)
elif all([os.path.isfile(f) for f in files2]):  ## foo.fa.ext
	for item in files2:
		path_to = os.path.join(new_dir, os.path.basename(item).replace(orig_prefix + file_extension, new_prefix + '.fasta'))
		os.symlink(item, path_to)
else:                                           ## create the index
	shell("samtools faidx {snakemake.params[0]}  {snakemake.input.fasta} > {snakemake.log}")
