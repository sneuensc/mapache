#!/usr/bin/env python

import os

pp = snakemake.params
print(pp)

if snakemake.params == "":
    os.system(f"cp {snakemake.input} {snakemake.output} 2>> {snakemake.log};")
else:
    os.system(
        f"bam trimBam {snakemake.input} {snakemake.output} {curparam} 2>> {snakemake.log};"
    )
