## indexing fasta
rule genome_index_bwa:
    """
    Indexing the genome for bwa
    """
    input:
        fasta="results/{folder}/{id_genome}.fasta"
    output:
        multiext("results/{folder}/{id_genome}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac")
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12)
    params:
        bwa_index_params = config.get("indexing", {}).get("bwa_params", ""),
        fasta=lambda wildcards: config.get("GENOME", {}).get(wildcards.id_genome, []),
    log:
        "results/logs/index/{folder}/bwa_index_{id_genome}.log"
    threads: 
        1
    conda:
    	"../envs/bwa.yaml"
    envmodules:
    	module_bwa
    message: "--- BWA INDEX  {input.fasta}"
    script:
    	"../scripts/bwa_indexing.py"
    
        
rule genome_index_bowtie2:
    """
    Indexing the genome for bowtie2
    """
    input:
        fasta="results/{folder}/{id_genome}.fasta"
    output:
        multiext("results/{folder}/{id_genome}.fasta", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12)
    params:
        bowtie2_index_params = config.get("indexing", {}).get("bowtie2_params", ""),
        fasta=lambda wildcards: config.get("GENOME", {}).get(wildcards.id_genome, []),
    log:
        "results/logs/index/{folder}/bowtie2_build_{id_genome}.log"
    threads: 
        1
    conda:
    	"../envs/bowtie2.yaml"
    envmodules:
    	module_bowtie2
    message: "--- BOWTIE2-BUILD  {input.fasta}"
    script:
        "../scripts/bowtie2_indexing.py"
                
rule samtools_index_fasta:
    """
    Indexing the genome with samtools faidx
    """
    input:
        fasta="results/{folder}/{id_genome}.fasta"
    output:
        "results/{folder}/{id_genome}.fasta.fai"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12)
    params:
        samtools_index_params = config.get("indexing", {}).get("samtools_params", ""),
        fasta=lambda wildcards: config.get("GENOME", {}).get(wildcards.id_genome, [])
    log:
        "results/logs/index/{folder}/samtools_fasta_{id_genome}.log"
    threads: 
    	1
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS_FAIDX {input.fasta}"
    script:
    	"../scripts/samtools_indexing.py"
        
rule genome_index_picard:
    """
    Indexing the genome with Picard CreateSequenceDictionary
    """
    input:
        fasta="results/{folder}/{id_genome}.fasta"
    output:
        "results/{folder}/{id_genome}.dict"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 12)
    params:
        picdard_index_params = config.get("indexing", {}).get("picard_params", ""),
        fasta=lambda wildcards: config.get("GENOME", {}).get(wildcards.id_genome, [])
    log:
        "results/logs/index/{folder}/picard_{id_genome}.log"
    threads: 
    	1
    conda:
    	"../envs/picard.yaml"
    envmodules:
    	module_picard
    message: "--- PICARD CreateSequenceDictionary {input.fasta}"
    script:
    	"../scripts/picard_indexing.py"
        


## indexing bam
rule samtools_index_bam:
    """
    Index bam file with samtools
    """
    input:
        "results/{folder}/{file}.bam"
    output:
        "results/{folder}/{file}.bai"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("indexing", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("indexing", attempt, 24)
    threads: 1
    log:
        "results/logs/{folder}/samtools_bam_{file}.log"
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	module_samtools
    message: "--- SAMTOOLS INDEX {input}"
    shell:
        "samtools index {input} {output} 2> {log};"