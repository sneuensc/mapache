#!/usr/bin/env python

import os

if snakemake.params.run and snakemake.params.number != 0:
    os.system(
        f"seqtk sample {snakemake.params.params} {snakemake.input} {snakemake.params.number} | gzip > {snakemake.output}"
    )
else:
    os.system(f"ln -srf {snakemake.input} {snakemake.output}")
