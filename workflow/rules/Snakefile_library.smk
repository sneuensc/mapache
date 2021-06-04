##########################################################################################
## all rules for libraries

def get_final_bam_library(id_sample, id_library, id_genome):
    if run_mapDamage_rescale:
        bam = "results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}.bam"
    elif extract_duplicates:
        bam = "results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}_mapped.bam"
    elif run_remove_duplicates:
        bam = "results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}.bam"
    else:
        bam = "results/02_library/00_merged_fastq/01_bam/{id_sample}/{id_library}.{id_genome}.bam"
    return (bam)


def get_mapDamage_bam(id_sample, id_library, id_genome, index = False):
	if extract_duplicates:
		bam="results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}_mapped.bam"
	elif run_remove_duplicates:
		bam="results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}.bam"
	else: 
		bam="results/02_library/00_merged_fastq/01_bam/{id_sample}/{id_library}.{id_genome}.bam"
	if index:
		bam = bam + ".bai"

	return (bam)


##########################################################################################
##########################################################################################

rule merge_bam_fastq2library:
    """
    Merge the bam files of the fastq step
    """
    input:
        lambda wildcards: ["results/01_fastq/04_final_fastq/01_bam/{SM}/{LB}/{FQ}.{genome}.bam".
        	format(FQ=FQ, LB=wildcards.id_library, SM=wildcards.id_sample, genome=wildcards.id_genome) 
            for FQ in db.loc[(db['LB']==wildcards.id_library) & (db['SM']==wildcards.id_sample)]['ID']]
    output:
        "results/02_library/00_merged_fastq/01_bam/{id_sample}/{id_library}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24)
    threads: 
        get_threads("merging", 4)
    log:
        "results/logs/02_library/00_merged_fastq/01_bam/{id_sample}/{id_library}.{id_genome}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS MERGE {input}"
    script:
    	"../scripts/merge_files.py"

        
rule merge_bam_fastq2library_low_qual:
    """
    Merge the bam files of the fastq step
    """
    input:
        lambda wildcards: ["results/01_fastq/04_final_fastq/01_bam_low_qual/{SM}/{LB}/{FQ}.{genome}.bam".
        	format(FQ=FQ, LB=wildcards.id_library, SM=wildcards.id_sample, genome=wildcards.id_genome) 
            for FQ in db.loc[(db['LB']==wildcards.id_library) & (db['SM']==wildcards.id_sample)]['ID']]
    output:
        "results/02_library/00_merged_fastq/01_bam_low_qual/{id_sample}/{id_library}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24)
    threads: 
    	get_threads("merging", 4)
    log:
        "results/logs/02_library/00_merged_fastq/01_bam_low_qual/{id_sample}/{id_library}.{id_genome}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS MERGE {input}"
    script:
    	"../scripts/merge_files.py"
        

rule remove_duplicates:
    """
    Remove duplicated mappings
    """
    input:
        "results/02_library/00_merged_fastq/01_bam/{id_sample}/{id_library}.{id_genome}.bam"
    output:
        bam="results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}.bam",
        stats="results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}.stats"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("markduplicates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("markduplicates", attempt, 24)
    params:
        rmdup_params = config.get("markduplicates", {}).get("params", "REMOVE_DUPLICATES=true"),
        PICARD = config.get("SOFTWARE", {}).get("picard_jar", "picard.jar")
    threads: 
    	get_threads("markduplicates", 4)
    log:
        "results/logs/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}.log"
    conda:
    	"../envs/picard.yaml"
    envmodules:
    	module_picard
    message: "--- MARKDUPLICATES {input}"
    shell:
    	"""
    	jar={params.PICARD}
        if [ "${{jar: -4}}" == ".jar" ]; then
        	java -XX:ParallelGCThreads={threads} -XX:+UseParallelGC -XX:-UsePerfData \
       			 -Xms{resources.memory}m -Xmx{resources.memory}m -jar {params.PICARD} MarkDuplicates \
        		INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.stats} {params.rmdup_params} \
        		ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> {log};
        else
        	{params.PICARD} MarkDuplicates \
        		INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.stats} {params.rmdup_params} \
        		ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> {log};
        fi
 		"""


rule samtools_extract_duplicates:
    """
    Extract duplicates of bam file
    """
    input:
        "results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}.bam"
    output:
        mapped="results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}_mapped.bam",
        dup="results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}_duplicates.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("markduplicates", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("markduplicates", attempt, 24)
    threads: 
    	get_threads("markduplicates", 4)
    log:
        "results/logs/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}_extract_duplicates.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS EXTRACT DUPLICATES {input}"
    shell:
        """
        samtools view -b --threads {threads} -F 1024 -U {output.dup} {input} > {output.mapped} 2> {log}
        """


rule mapDamage_stats:
    """
    Run mapDamage to quantify the deamination pattern
    """
    input:
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        bam=lambda wildcards: get_mapDamage_bam(wildcards.id_sample, wildcards.id_library, wildcards.id_genome)
    output:
        #directory("results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}_results_mapDamage"),
        report("results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}_results_mapDamage/Fragmisincorporation_plot.pdf", category="Damage pattern"),
        report("results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}_results_mapDamage/Length_plot.pdf", category="Read length distribution")
    log:
        "results/logs/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}_stats.log"
    threads: 1
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24)
    params:
        prefix="results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}",
        mapDamage_params = config.get("mapdamage", {}).get("params", "") 
    conda:
    	"../envs/mapdamage.yaml"
    message: "--- MAPDAMAGE {input.bam}"
    shell:
        """
        mapDamage -i {input.bam} -r {input.ref} -d {params.prefix}_results_mapDamage \
        {params.mapDamage_params} --merge-reference-sequences 2> {log};
        """


rule mapDamage_rescale:
    """
    Run mapDamage to rescale bam file
    """
    input:
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        bam=lambda wildcards: get_mapDamage_bam(wildcards.id_sample, wildcards.id_library, wildcards.genome),
        stats="results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}_results_mapDamage/Fragmisincorporation_plot.pdf"
    output:
        "results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapdamage", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapdamage", attempt, 24)
    log:
        "results/logs/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}_rescale.log"
    threads: 1
    params:
        prefix="results/02_library/02_rescaled/01_mapDamage/{id_sample}/{id_library}.{id_genome}",
        mapDamage_params = config.get("mapdamage", {}).get("params", "")
    conda:
    	"../envs/mapdamage.yaml"
    message: "--- MAPDAMAGE {input.bam}"
    shell:
    	"""
    	mapDamage -i {input.bam} -r {input.ref} -d {params.prefix}_results_mapDamage \
        {params.mapDamage_params} --merge-reference-sequences --rescale-only --rescale-out {output} 2>> {log};
        """
        

##########################################################################################
rule get_final_library:
	"""
	Get the final bam file of the library part
	"""
	input:
		lambda wildcards: get_final_bam_library(wildcards.id_sample, wildcards.id_library, wildcards.id_genome)
	output:
		"results/02_library/03_final_library/01_bam/{id_sample}/{id_library}.{id_genome}.bam"
	message: "--- GET FINAL BAM {input} (LIBRARY LEVEL)"
	run:
		symlink_rev(input, output)
	
rule get_final_library_low_qual:
	"""
	Get the final bam file of the library part
	"""
	input:
		"results/02_library/00_merged_fastq/01_bam_low_qual/{id_sample}/{id_library}.{id_genome}.bam"
	output:
		"results/02_library/03_final_library/01_bam_low_qual/{id_sample}/{id_library}.{id_genome}.bam"	
	message: "--- GET FINAL LOW_QUAL BAM {input} (LIBRARY LEVEL)"
	run:
		symlink_rev(input, output)
		
rule get_final_library_lduplicate:
	"""
	Get the final duplicate bam file of the library part
	"""
	input:
		"results/02_library/01_duplicated/01_rmdup/{id_sample}/{id_library}.{id_genome}_duplicates.bam"
	output:
		"results/02_library/03_final_library/01_bam_duplicate/{id_sample}/{id_library}.{id_genome}.bam"	
	message: "--- GET FINAL DUPLICATE BAM {input} (LIBRARY LEVEL)"
	run:
		symlink_rev(input, output)

