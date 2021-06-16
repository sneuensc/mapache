##########################################################################################
## all rules for the samples
##########################################################################################
localrules: get_final_bam ## executed locally on a cluster
ruleorder:  merge_bam_library2sample_low_qual  > get_final_bam 

rule merge_bam_library2sample:
    """
    Merge the bam files of the library step
    """
    input:
    	mapped=lambda wildcards: [("results/02_library/03_final_library/01_bam/{SM}/{LB}.{genome}.bam").
    		format(LB=LB, SM=wildcards.id_sample, genome=wildcards.id_genome) 
    		for LB in samples[wildcards.id_sample].keys()]
    output:
        "results/03_sample/00_merged_library/01_bam/{id_sample}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24)
    threads: 
    	get_threads("merging", 4)
    log:
        "results/logs/03_sample/00_merged_library/01_bam/{id_sample}.{id_genome}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS MERGE merge_bam_library2sample {input}"
    script:
    	"../scripts/merge_files.py"


rule merge_bam_library2sample_low_qual:
    """
    Merge the low quality bam files of the library step
    """
    input:
    	lambda wildcards: [("results/02_library/03_final_library/01_bam_low_qual/{SM}/{LB}.{genome}.bam").
    		format(LB=LB, SM=wildcards.id_sample, genome=wildcards.id_genome) 
    		for LB in samples[wildcards.id_sample].keys()]
    output:
        "results/03_sample/00_merged_library/01_bam_low_qual/{id_sample}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24)
    threads: 
    	get_threads("merging", 4)
    log:
        "results/logs/03_sample/00_merged_library/01_bam_low_qual/{id_sample}.{id_genome}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS MERGE merge_bam_library2sample_low_qual {input}"
    script:
    	"../scripts/merge_files.py"


rule merge_bam_library2sample_duplicates:
    """
    Merge the duplicate bam files of the library step
    """
    input:
        low_qual=lambda wildcards: get_bams_of_sample_low_cov("library_rmdup", "_duplicates.bam", wildcards.id_sample, wildcards.id_genome)
    output:
        "results/03_sample/00_merged_library/01_bam_duplicate/{id_sample}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("merging", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("merging", attempt, 24)
    threads: 
    	get_threads("merging", 4)
    log:
        "results/logs/03_sample/00_merged_library/01_bam_duplicate/{id_sample}.{id_genome}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS MERGE merge_bam_library2sample_duplicates {input}"
    script:
    	"../scripts/merge_files.py"


rule realign:
    """
    Realign sequence around indels.
    """
    input:
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        fai="results/00_reference/{id_genome}/{id_genome}.fasta.fai",
        dict="results/00_reference/{id_genome}/{id_genome}.dict",
        bam="results/03_sample/00_merged_library/01_bam/{id_sample}.{id_genome}.bam",
        bai="results/03_sample/00_merged_library/01_bam/{id_sample}.{id_genome}.bai"
    output:
        bam="results/03_sample/01_realigned/01_realign/{id_sample}.{id_genome}.bam",
        intervals="results/03_sample/01_realigned/01_realign/{id_sample}.{id_genome}.intervals",
        bai="results/03_sample/01_realigned/01_realign/{id_sample}.{id_genome}.bai"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("realign", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("realign", attempt, 24)
    threads: 
    	get_threads("realign", 4)
    params:
    	GATK = config.get("SOFTWARE", {}).get("gatk3_jar", "GenomeAnalysisTK.jar")  
    log:
        "results/logs/03_sample/01_realigned/01_realign/{id_sample}.{id_genome}.log"
    conda:
    	"../envs/gatk3.yaml"
    envmodules:
    	module_gatk3	
    message: "--- GATK INDELREALIGNER {input.bam}"
    shell:
        """
     	jar={params.GATK};
        if [ "${{jar: -4}}" == ".jar" ]; then
       		java -Djava.io.tmpdir=/tmp/ -XX:ParallelGCThreads={threads} -XX:+UseParallelGC \
        		-XX:-UsePerfData -Xms15000m -Xmx15000m -jar {params.GATK} \
        		-I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.intervals} 2> {log}; \
        	java -Djava.io.tmpdir=/tmp/ -XX:ParallelGCThreads={threads} -XX:+UseParallelGC \
        		-XX:-UsePerfData -Xms15000m -Xmx15000m -jar {params.GATK} \
        		-I {input.bam} -T IndelRealigner -R {input.ref} -targetIntervals \
        		{output.intervals} -o {output.bam} 2>> {log};
        else
       		{params.GATK} -I {input.bam} -R {input.ref} -T RealignerTargetCreator -o {output.intervals} 2> {log}; \
        	{params.GATK} -I {input.bam} -T IndelRealigner -R {input.ref} -targetIntervals \
        		{output.intervals} -o {output.bam} 2>> {log};
        fi
        """

            
rule samtools_calmd:
    """
    Recompute the md flag.
    """
    input:
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        bam=lambda wildcards: get_md_flag_bam(wildcards.id_sample, wildcards.id_genome)
    output:
        "results/03_sample/02_md_flag/01_md_flag/{id_sample}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("calmd", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("calmd", attempt, 24)
    threads: 
    	get_threads("calmd", 4)
    log:
        "results/logs/03_sample/02_md_flag/01_md_flag/{id_sample}.{id_genome}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS CALMD {input.bam}"
    shell:
        """
        samtools calmd --threads {threads} {input.bam} {input.ref} 2> {log} | samtools view -bS - > {output}
        """
        
        
rule get_final_bam:
    """
    Get the final bam files 
    """
    input:
        lambda wildcards: get_final_bam(wildcards.id_sample, wildcards.id_genome)
    output:
        "results/03_sample/03_final_sample/01_bam/{id_sample}.{id_genome}.bam"
    threads: 1
    message: "--- SIMLINKK FINAL BAM"
    run:
        symlink_rev(input, output)
        
rule get_final_bam_low_qual:
    """
    Get the final bam files 
    """
    input:
        "results/03_sample/00_merged_library/01_bam_low_qual/{id_sample}.{id_genome}.bam"
    output:
        "results/03_sample/03_final_sample/01_bam_low_qual/{id_sample}.{id_genome}.bam"
    threads: 1
    message: "--- SIMLINKK FINAL LOW_QUAL BAM"
    run:
        symlink_rev(input, output)
        
rule move_final_bam_duplicate:
    """
    Get the final bam files
    """
    input:
        "results/03_sample/00_merged_library/01_bam_duplicate/{id_sample}.{id_genome}.bam"
    output:
        "results/03_sample/03_final_sample/01_bam_duplicate/{id_sample}.{id_genome}.bam"
    threads: 1
    message: "--- SIMLINKK FINAL DUPLICATE BAM"
    run:
        symlink_rev(input, output)
        

        
##########################################################################################
