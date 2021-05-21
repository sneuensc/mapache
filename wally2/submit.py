#!/usr/bin/env python3

'''
Submission script for Snakemake.
The entire submission process is handled by the SlurmScheduler object.
'''

from scheduler import SlurmScheduler

if __name__ == '__main__':
    scheduler = SlurmScheduler()
    scheduler.submit()
