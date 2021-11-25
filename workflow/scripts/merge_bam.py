#!/usr/bin/env python

import os

## Merge bam files. If it is a single file, just symlink the file.
if len(snakemake.input)>1:
	os.system(f'samtools merge -f --threads {snakemake.threads} {snakemake.output} {snakemake.input} 2> {snakemake.log};')
else: 
	os.system(f'ln -srf {snakemake.input} {snakemake.output}')
