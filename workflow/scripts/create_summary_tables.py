#!/usr/bin/env python

import os
import pandas as pd
import numpy as np

## input files
fastqc_orig = snakemake.input.fastqc_orig[0]
fastqc_trim = snakemake.input.fastqc_trim[0]
flagstat_fastq_final = snakemake.input.flagstat_fastq_final
flagstat_library_final = snakemake.input.flagstat_library_final

## output files
fastq_stats = snakemake.output.fastq_stats
library_stats = snakemake.output.library_stats
sample_stats = snakemake.output.sample_stats    

## params
delim = snakemake.params.delim   
sample_file = snakemake.params.sample_file

## log
log = snakemake.log


def read_flagstat(file, extract, name):
	d1 = pd.read_csv(file, sep="\t").rename({extract: name}, axis=1)[['Sample', name]] 
	d1[name] = d1[name].astype('int64')
	d2 = pd.DataFrame(d1['Sample'].str.split('|').tolist(), columns = ['SM','LB','type','ID'])
	d2['SM'] = d2['SM'].str.strip()
	d2['LB'] = d2['LB'].str.strip()
	d2['ID'] = d2['ID'].str.strip()   
	if "_R1_data" == os.path.dirname(file)[-8:]:
		d2['ID'] = d2['ID'].map(lambda x: str(x)[:-3])
	d2['ID'] = d2['ID'].map(lambda x: x.rstrip(".{}_flagstat".format(snakemake.wildcards.id_genome)))
	dd = pd.concat([d2[['SM','LB','ID']], d1[name]], axis=1)
	return (dd)


## sample file    
db_fastq = pd.read_csv(sample_file, sep=delim)[['SM', 'LB' ,'ID']]\
.sort_values(['SM', 'LB' ,'ID'], axis=0, ascending=True)\
.reset_index(drop=True)

## fastq file
db_orig = read_flagstat(fastqc_orig, 'Total Sequences', 'reads_raw')
db_trim = read_flagstat(fastqc_trim, 'Total Sequences', 'reads_trim')
db_map = read_flagstat(flagstat_fastq_final, 'mapped_passed', 'mapping')

db_fastq = pd.merge(db_fastq, db_orig, how="left", on = ['SM','LB','ID'])
db_fastq = pd.merge(db_fastq, db_trim, how="left", on = ['SM','LB','ID'])
db_fastq = pd.merge(db_fastq, db_map, how="left", on = ['SM','LB','ID'])
db_fastq['trim_prop'] = np.where(db_fastq['reads_raw'].ne(0), db_fastq['reads_trim'].astype(float) / db_fastq['reads_raw'], np.NaN) 
db_fastq['endo_prop'] = np.where(db_fastq['reads_raw'].ne(0), db_fastq['mapping'].astype(float) / db_fastq['reads_raw'], np.NaN) 

db_fastq[['ID', 'LB', 'SM', 'reads_raw', 'reads_trim', 'trim_prop', 'mapping', 'endo_prop']]\
.round({'reads_trim': 0, 'trim_prop': 4, 'mapping': 0, 'endo_prop': 4})\
.to_csv(fastq_stats, index = None, header=True)

## library stats
db_library = db_fastq.groupby(['SM','LB'])[['reads_raw', 'reads_trim', 'mapping']].apply(lambda x: x.sum(skipna=False)).reset_index()
db_lib = read_flagstat(flagstat_library_final, 'mapped_passed', 'mapping_final').drop('ID', 1)

db_library = pd.merge(db_library, db_lib, how="left", on = ['SM','LB'])
db_library['duplicates'] = db_library['mapping'] - db_library['mapping_final']
db_library['trim_prop'] = np.where(db_library['reads_raw'].ne(0), db_library['reads_trim'].astype(float) / db_library['reads_raw'], np.NaN) 
db_library['endo_prop'] = np.where(db_library['reads_raw'].ne(0), db_library['mapping'].astype(float) / db_library['reads_raw'], np.NaN) 
db_library['endo_final_prop'] = np.where(db_library['reads_raw'].ne(0), db_library['mapping_final'].astype(float) / db_library['reads_raw'], np.NaN) 
db_library['duplicates_prop'] = np.where(db_library['mapping'].ne(0), db_library["duplicates"].astype(float) / db_library["mapping"], np.NaN) 

db_library[['LB', 'SM','reads_raw', 'reads_trim', 'trim_prop', 'mapping', 'endo_prop', 'duplicates', 'duplicates_prop', 'mapping_final', 'endo_final_prop']]\
.round({'reads_raw': 0, 'reads_trim': 0, 'trim_prop': 4, 'mapping': 0, 'endo_prop': 4, 'duplicates': 0, 'duplicates_prop': 4, 'mapping_final': 0, 'endo_final_prop': 4})\
.to_csv(library_stats, index = None, header=True)

## sample stats
db_sample = db_library.groupby('SM')[['reads_raw', 'reads_trim', 'mapping', 'mapping_final', 'duplicates']].apply(lambda x: x.sum(skipna=False)).reset_index()
db_sample['trim_prop'] = np.where(db_sample['reads_raw'].ne(0), db_sample['reads_trim'].astype(float) / db_sample['reads_raw'], np.NaN) 
db_sample['endo_prop'] = np.where(db_sample['reads_raw'].ne(0), db_sample['mapping'].astype(float) / db_sample['reads_raw'], np.NaN) 
db_sample['endo_final_prop'] = np.where(db_sample['reads_raw'].ne(0), db_sample['mapping_final'].astype(float) / db_sample['reads_raw'], np.NaN) 
db_sample['duplicates_prop'] = np.where(db_sample['mapping'].ne(0), db_sample["duplicates"].astype(float) / db_sample["mapping"], np.NaN) 

db_sample[['SM','reads_raw', 'reads_trim', 'trim_prop', 'mapping', 'endo_prop', 'duplicates', 'duplicates_prop', 'mapping_final', 'endo_final_prop']]\
.round({'reads_raw': 0, 'reads_trim': 0, 'trim_prop': 4, 'mapping': 0, 'endo_prop': 4, 'duplicates': 0, 'duplicates_prop': 4, 'mapping_final': 0, 'endo_final_prop': 4})\
.to_csv(sample_stats, index = None, header=True)
