#!/usr/bin/env python

import os
	
if len(snakemake.input)>1:
	os.system(f'samtools merge -f --threads {snakemake.threads} {snakemake.output} {snakemake.input} 2> {snakemake.log};')
else: 
	os.system(f'ln -srf {snakemake.input} {snakemake.output}')


#        """
# 		nb=$(echo {input} | wc -w);
#     	if [ $"nb" -eq 1 ]
#     	then
# 			ln -srf {input} {output};
# 		else
# 			samtools merge -f --threads {threads} {output} {input} 2> {log};			
# 		fi
#        """ 	