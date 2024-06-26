#!/usr/bin/env python

## This script tests if a given index is present next to the original file.
## If present the index is symlinked, otherwise created.

import subprocess
import os

## get the original directory
orig_fasta = snakemake.input.orig
orig_prefix = os.path.splitext(orig_fasta)[0]

## get the extensions
out_files = list(snakemake.output)
in_fasta = snakemake.input.fasta
ext = [s.replace(in_fasta, "") for s in out_files]

## check if the indexes are present
files1 = [orig_fasta + s for s in ext]  ## foo.fasta.ext
files2 = [orig_prefix + s for s in ext]  ## foo.ext
if all([os.path.isfile(f) for f in files1]):  ## foo.fasta.ext
    for idx, item in enumerate(files1):
        os.symlink(item, out_files[idx])
elif all([os.path.isfile(f) for f in files2]):  ## foo.ext
    for idx, item in enumerate(files2):
        os.symlink(item, out_files[idx])
else:  ## index is not present: create the index
    cmd = eval(snakemake.params.cmd)
    # print(cmd)
    subprocess.run(cmd, shell=True)
