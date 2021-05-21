#!/usr/bin/env python3
"""
Snakemake SLURM submit script.
"""
from snakemake.utils import read_job_properties

import slurm_utils
from CookieCutter import CookieCutter

# cookiecutter arguments
SBATCH_DEFAULTS = CookieCutter.SBATCH_DEFAULTS
CLUSTER = CookieCutter.get_cluster_option()
CLUSTER_CONFIG = CookieCutter.CLUSTER_CONFIG
ADVANCED_ARGUMENT_CONVERSION = CookieCutter.get_advanced_argument_conversion()

RESOURCE_MAPPING = {
    "time": ("time", "runtime", "walltime"),
    "mem": ("mem", "mem_mb", "ram", "memory"),
    "mem-per-cpu": ("mem-per-cpu", "mem_per_cpu", "mem_per_thread"),
    "nodes": ("nodes", "nnodes"),
}

# parse job
jobscript = slurm_utils.parse_jobscript()
job_properties = read_job_properties(jobscript)

sbatch_options = {}
cluster_config = slurm_utils.load_cluster_config(CLUSTER_CONFIG)

# 1) sbatch default arguments and cluster
sbatch_options.update(slurm_utils.parse_sbatch_defaults(SBATCH_DEFAULTS))
sbatch_options.update(slurm_utils.parse_sbatch_defaults(CLUSTER))

# 2) cluster_config defaults
sbatch_options.update(cluster_config["__default__"])

# 3) Convert resources (no unit conversion!) and threads
sbatch_options.update(
    slurm_utils.convert_job_properties(job_properties, RESOURCE_MAPPING)
)

# 4) cluster_config for particular rule
sbatch_options.update(cluster_config.get(job_properties.get("rule"), {}))

# 5) cluster_config options
sbatch_options.update(job_properties.get("cluster", {}))

# 6) Advanced conversion of parameters
if ADVANCED_ARGUMENT_CONVERSION:
    sbatch_options = slurm_utils.advanced_argument_conversion(sbatch_options)

# 7) Format pattern in snakemake style
sbatch_options = slurm_utils.format_values(sbatch_options, job_properties)

# ensure sbatch output dirs exist
for o in ("output", "error"):
    slurm_utils.ensure_dirs_exist(sbatch_options[o]) if o in sbatch_options else None

# submit job and echo id back to Snakemake (must be the only stdout)
print(slurm_utils.submit_job(jobscript, **sbatch_options))
