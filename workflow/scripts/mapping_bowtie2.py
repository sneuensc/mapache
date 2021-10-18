#!/usr/bin/env python

import os

if len(input.fastq) == 2:
	os.system(f'bowtie2 {snakemake.params.bowtie2_params} -p {snakemake.threads} -x {snakemake.input.ref} \
		  -1 {snakemake.input.fastq[0]} -2 {snakemake.input.fastq[1]} 2> {snakemake.log} | \
		  samtools view -bS - > {snakemake.output}')
else:
	os.system('fbowtie2 {snakemake.params.bowtie2_params} -p {snakemake.threads} -x {snakemake.input.ref} \
		  -U {snakemake.input.fastq} 2> {snakemake.log} | samtools view -bS - > {snakemake.output}')
