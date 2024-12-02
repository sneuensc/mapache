#!/usr/bin/env python

import subprocess
import os

## change sample name in RG
def reheader_bam(in_file, out_file, sm):
    header = f"{out_file}_header.sam"
    cmd = f"samtools view -H {in_file} | sed 's/<SM:[[:space:]]*[^[:space:]]*>/SM:{sm}/g' > {header}; \
            samtools reheader -P {header} {in_file} > {out_file};"
    subprocess.run(cmd, shell=True)
    os.remove(header)

out_file = snakemake.output[0]

## multiple bam files: merge them
if len(list(snakemake.input)) > 1:
    cmd = f"samtools merge -f --threads {snakemake.threads} {out_file} {snakemake.input} 2> {snakemake.log};"
    subprocess.run(cmd, shell=True)

    ## rename header if needed
    if snakemake.params.sm_changed:
        temp_file = f"{out_file}.bam"
        reheader_bam(out_file, temp_file, snakemake.wildcards.sm)
        os.rename(temp_file, out_file)

## a single bam file: just overwrite read group
elif snakemake.params.sm_changed:
    reheader_bam(snakemake.input[0], out_file, snakemake.wildcards.sm)

else:
    LOGGER.error(
        f"ERROR: Script merge_bam.py should never pass here!"
    )
    os._exit(1)