#!/usr/bin/env python

import os
import pandas as pd

## concatenate depths
df = pd.DataFrame(columns=['Sample','Library','AvgReadDepth_tot', 
						   'AvgReadDepth_MT', 'Nseqs', 'NchrY+NchrX', 
						   'NchrY', 'R_y', 'SE', '95%CI', 'Sex', 
						   'AvgReadLength'])

for file in snakemake.input:
	## read the file
	to_store = []
	for cnt, row in enumerate(open(file)):
		to_store.append(row)
	
	if len(to_store) == 0:
		continue
	
	## extract info
	# Sample id, genome coverage, mt coverage
	a = [os.path.basename(file).replace('_depth.txt', '').replace(".{}".format(wildcards.id_genome), '')]
	## AvgReadDepth_tot
	curLine = list(filter(lambda x:'Total' in x, to_store))
	a.append(curLine[0].split()[1] if len(curLine)==1 else 'NA')
	## AvgReadDepth_MT
	curLine = list(filter(lambda x:'M' in x, to_store))
	a.append(curLine[0].split()[1] if len(curLine)==1 else 'NA')
	## 'Nseqs', 'NchrY+NchrX', 'NchrY', 'R_y', 'SE', '95%CI', 'Sex'
	curLine = list(filter(lambda x:'Nseqs' in x, to_store))
	lineNb = to_store.index(curLine[0]) if len(curLine)==1 else 0
	if lineNb is not 0:
		a += to_store[lineNb+1].replace('\n', '').split('\t')
	else:
		a += ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']
	
	##     AvgReadLength
	curLine = list(filter(lambda x:'Average length of mapped reads' in x, to_store))
	a.append([curLine[0].split()[-1]][0] if len(curLine)==1 else 'NA')
	df.loc[len(df)] = a    
	
df.to_csv(snakemake.output.library_depth, index = None, header=True)
