##########################################################################################
## all rules for fastq files
##########################################################################################
## get/rename reference and fastq files

localrules: get_fastq, get_fasta ## executed locally on a cluster
#ruleorder: genome_index_bwa > get_fasta

## all rules for fastq files
rule get_fastq:
    """
    Symlink and rename all fastq files to a common folder (makes the DAG readable)
    """
    input:
        get_fastq_of_ID
    output:
        "{folder}/00_reads/01_files_orig/{id_sample}/{id_library}/{id_fastq}.fastq.gz"
    threads: 1
    message: "--- GET FASTQ FILES  {input}"            
    shell:
    	"""
    	ln -srf {input} {output}
    	"""


## all rules for fastq files
rule get_fasta:
    """
    Symlink and rename the reference (.fasta/.fa) to a new folder.
    """
    input:
       lambda wildcards: get_param3('genome', wildcards.id_genome, 'fasta', '')
    output:
        "results/00_reference/{id_genome}/{id_genome}.fasta",
    threads: 1
    message: "--- GET REFERENCE  {input}" 
    shell:
    	"""
    	ln -srf {input} {output}
    	"""



##########################################################################################
## trimming

ruleorder: adapter_removal_pe > adapter_removal_se


rule adapter_removal_se:
    """
    Remove adapter and low quality bases at the edges
    """
    input:
        "{folder}/00_reads/01_files_orig/{id_sample}/{id_library}/{id_fastq}.fastq.gz"
    output:
        fastq = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.fastq.gz",
        discard = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.discarded.gz",
        setting = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.settings"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("adapterremoval", attempt, 4),
        runtime = lambda wildcards, attempt: get_runtime_alloc("adapterremoval", attempt, 24)
    params:
        get_param2("adapterremoval", "params", "--minlength 30 --trimns --trimqualities")
    log:
        "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.log"
    threads: 
    	get_threads("adapterremoval", 4)
    conda:
    	"../envs/adapterremoval.yaml"
    message: "--- ADAPTERREMOVAL  {input}"            
    shell:
        """
        out={output.fastq};
        AdapterRemoval --threads {threads} {params} --file1 {input} \
                --basename ${{out%%.fastq.gz}} --gzip \
                --output1 {output.fastq} 2> {log};
        """


rule adapter_removal_pe:
    """
    Remove adapter and low quality bases at the edges
    """
    input:
        R1 = "{folder}/00_reads/01_files_orig/{id_sample}/{id_library}/{id_fastq}_R1.fastq.gz",
        R2 = "{folder}/00_reads/01_files_orig/{id_sample}/{id_library}/{id_fastq}_R2.fastq.gz"
    output:
        R1 = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}_R1.fastq.gz",
        R2 = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}_R2.fastq.gz",
        singleton = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.singleton.truncated.gz",
        discarded = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.discarded.gz",
        settings = "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.settings"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("adapterremoval", attempt, 4),
        runtime = lambda wildcards, attempt: get_runtime_alloc("adapterremoval", attempt, 24)
    params:
        get_param2("adapterremoval", "params", "--minlength 30 --trimns --trimqualities")
    log:
        "{folder}/01_trimmed/01_files_trim/{id_sample}/{id_library}/{id_fastq}.log"
    threads: 
    	get_threads("adapterremoval", 4)
    conda:
    	"../envs/adapterremoval.yaml"
    message: "--- ADAPTERREMOVAL {input.R1} {input.R2}"            
    shell:
        """
        out={output.R1};
        AdapterRemoval --threads {threads} {params} --file1 {input.R1} \
                --file2 {input.R2} --basename ${{out%%_R1.fastq.gz}} --gzip \
                --output1 {output.R1} --output2 {output.R2} 2> {log};
        """


rule adapter_removal_collapse:
    """
    Remove adapter and low quality bases at the edges and collapse paired-end reads
    """
    input:
        R1 = "{folder}/00_reads/01_files_orig/{id_sample}/{id_library}/{id_fastq}_R1.fastq.gz",
        R2 = "{folder}/00_reads/01_files_orig/{id_sample}/{id_library}/{id_fastq}_R2.fastq.gz"
    output:
        R = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}.fastq.gz",
        trunc = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}_truncated.fastq.gz",
        R1 = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}_R1.fastq.gz",
        R2 = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}_R2.fastq.gz",
        strunc = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}.singleton.truncated.gz",
        disc = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}.discarded.gz",
        settingd = "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}.settings"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("adapterremoval", attempt, 4),
        runtime = lambda wildcards, attempt: get_runtime_alloc("adapterremoval", attempt, 24)
    params:
        get_param2("adapterremoval", "params", "--minlength 30 --trimns --trimqualities")
    log:
        "{folder}/01_trimmed/01_files_trim_collapsed/{id_sample}/{id_library}/{id_fastq}.log"
    threads: 
    	get_threads("adapterremoval", 4)
    conda:
    	"../envs/adapterremoval.yaml"
    message: "--- ADAPTERREMOVAL {input.R1} {input.R2}"            
    shell:
        """
        out={output.R};
        AdapterRemoval --threads {threads} {params} --file1 {input.R1} \
                --file2 {input.R2} --basename ${{out%%.fastq.gz}} --gzip \
                --output1 {output.R1} --output2 {output.R2} \
                --outputcollapsed {output.R} \
                --outputcollapsedtruncated {output.trunc} 2> {log};
        """



##########################################################################################
## mapping
ruleorder: mapping_bwa_aln_pe > mapping_bwa_aln_se

rule mapping_bwa_aln_se:
    """
    Align reads to the reference
    """
    input:
        multiext("results/00_reference/{id_genome}/{id_genome}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac"),
        ref = "results/00_reference/{id_genome}/{id_genome}.fasta",
        fastq = get_fastq_for_mapping
    output:
        "{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}.sai"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime = lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24)
    params:
        config.get("mapping", {}).get("bwa_aln_params", "-l 1024")
    log:
        "{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
    threads: 
    	get_threads("mapping", 4)
    conda:
    	"../envs/bwa.yaml"
    envmodules:
    	module_bwa
    message: "--- BWA ALN  {input.fastq}"
    shell:
        """
        bwa aln {params} -t {threads} {input.ref} -f {output} {input.fastq} 2> {log}
        """

rule mapping_bwa_aln_pe:
    """
    Align reads to the reference
    """
    input:
        multiext("results/00_reference/{id_genome}/{id_genome}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac"),
        ref = "results/00_reference/{id_genome}/{id_genome}.fasta",
        fastq = get_fastq_for_mapping_pe
    output:
        "{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}_R{id_read}.sai"
    resources:
        memory = lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime = lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24)
    params:
        config.get("mapping", {}).get("bwa_aln_params", "-l 1024") 
    log:
        "{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}_R{id_read}.log"
    threads: 
    	get_threads("mapping", 4)
    conda:
    	"../envs/bwa.yaml"
    envmodules:
    	module_bwa
    message: "--- BWA ALN  {input.fastq}"
    shell:
        """
        bwa aln {params} -t {threads} {input.ref} -f {output} {input.fastq} 2> {log}
        """


rule mapping_bwa_samse:
    """
    Creates bam file from sai file for SE reads
    """
    input:
        multiext("results/00_reference/{id_genome}/{id_genome}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac"),
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        fastq=get_fastq_for_mapping,
        sai="{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}.sai"
    output:
        "{folder}/02_mapped/02_bwa_samse/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24)
    params:
        ID="{id_fastq}",
        LB=lambda wildcards: get_from_sample_file("LB", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        SM=lambda wildcards: get_from_sample_file("SM", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        PL=lambda wildcards: get_from_sample_file("PL", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        bwa_samse_params = config.get("mapping", {}).get("bwa_samse_params", "-n 3")
    log:
        "{folder}/02_mapped/02_bwa_samse/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
    threads: 1
    conda:
    	"../envs/bwa.yaml"
    envmodules:
    	"gcc bwa/0.7.17",
    	"gcc samtools/1.4"
    message: "--- BWA SAMSE  {input.fastq}"
    shell:
        """
        (bwa samse {params.bwa_samse_params} \
         -r \"@RG\\tID:{params.ID}\\tLB:{params.LB}\\tSM:{params.SM}\\tPL:{params.PL}\" \
         {input.ref} {input.sai} {input.fastq} | samtools view -Sb > {output}) 2> {log}
         """


rule mapping_bwa_sampe:
    """
    Creates bam file from sai file for PE reads
    """
    input:
        multiext("results/00_reference/{id_genome}/{id_genome}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac"),
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        fastq=get_fastq_for_mapping,	## should get both pairs
        sai1="{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}_R1.sai",
        sai2="{folder}/02_mapped/01_bwa_aln/{id_sample}/{id_library}/{id_fastq}.{id_genome}_R2.sai"
    output:
        "{folder}/02_mapped/02_bwa_sampe/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24)
    params:
        ID="{id_fastq}",
        LB=lambda wildcards: get_from_sample_file("LB", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        SM=lambda wildcards: get_from_sample_file("SM", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        PL=lambda wildcards: get_from_sample_file("PL", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        bwa_samse_params = config.get("mapping", {}).get("bwa_samse_params", "-n 3")
    log:
        "{folder}/02_mapped/02_bwa_sampe/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
    threads: 1
    conda:
    	"../envs/bwa.yaml"
    envmodules:
    	module_bwa,
    	module_samtools
    message: "--- BWA SAMPE {input.fastq}"
    shell:
        """
        (bwa sampe {params.bwa_samse_params} \"
         	-r \"@RG\\tID:{params.ID}\\tLB:{params.LB}\\tSM:{params.SM}\\tPL:{params.PL}\" \
         	{input.ref} {input.sai1} {input.sai2} {input.fastq} | \
         	samtools view -Sb > {output}) 2> {log}
         """


rule mapping_bwa_mem:
    """
    Map reads to genome using bwa mem
    """
    input:
        multiext("results/00_reference/{id_genome}/{id_genome}.fasta", ".sa", ".amb", ".ann", ".bwt", ".pac"),
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        fastq=get_fastq_for_mapping
    output:
        "{folder}/02_mapped/02_bwa_mem/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24)
    params:
        ID="{id_fastq}",
        LB=lambda wildcards: get_from_sample_file("LB", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        SM=lambda wildcards: get_from_sample_file("SM", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        PL=lambda wildcards: get_from_sample_file("PL", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        bwa_mem_params = config.get("mapping", {}).get("bwa_mem_params", "")
    log:
        "{folder}/02_mapped/02_bwa_mem/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
    threads: 
    	get_threads("mapping", 4)
    conda:
    	"../envs/bwa.yaml"
    envmodules:
    	module_bwa
    message: "--- BWA MEM {input.fastq}"
    shell:
        """
        bwa mem {params.bwa_mem_params} -t {threads} \
    		-R \"@RG\\tID:{params.ID}\\tLB:{params.LB}\\tSM:{params.SM}\\tPL:{params.PL}\" \
        	{input.ref} ${input.fastq} > {output} 2> {log};
        """


rule mapping_bowtie2:
    """
    Map reads to genome using bowtie2
    """
    input:
        multiext("results/00_reference/{id_genome}/{id_genome}.fasta", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        ref="results/00_reference/{id_genome}/{id_genome}.fasta",
        fastq=get_fastq_for_mapping
    output:
        "{folder}/02_mapped/02_bwa_bowtie2/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mapping", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("mapping", attempt, 24)
    params:
        ID="{id_fastq}",
        LB=lambda wildcards: get_from_sample_file("LB", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        SM=lambda wildcards: get_from_sample_file("SM", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        PL=lambda wildcards: get_from_sample_file("PL", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq)[0],
        bowtie2_params = config.get("mapping", {}).get("bowtie2_params", "")
    log:
        "{folder}/02_mapped/02_bwa_bowtie2/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
    threads: 
    	get_threads("mapping", 4)
    conda:
    	"../envs/bowtie2.yaml"
    envmodules:
    	module_bowtie2,
    	module_samtools
    message: "--- BOWTIE2 {input.fastq}"
    script:
    	"../scripts/mapping_bowtie2.py"


##########################################################################################
## sorting

rule samtools_sort:
    """
    Sort bam file with samtools
    """
    input:
        get_bam_for_sorting
    output:
        "{folder}/02_mapped/03_bam_sort/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("sorting", attempt, 4),
        runtime=lambda wildcards, attempt: get_runtime_alloc("sorting", attempt, 24)
    log:
        "{folder}/02_mapped/03_bam_sort/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
    threads: 
    	get_threads("sorting", 4)
    conda:
    	"../envs/samtools.yaml"
    envmodules:
    	"gcc samtools/1.4"
    message: "--- SAMTOOLS SORT {input}"
    shell:
        """
        samtools sort --threads {threads} {input} > {output} 2> {log}
        """
        

##########################################################################################
## filtering

if save_low_qual:
	rule samtools_filter:
		"""
		Filter mappings following quality and keeping the low quality mappings
		"""
		input:
			"{folder}/02_mapped/03_bam_sort/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
		output:
			mapped="{folder}/03_filtered/01_bam_filter/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam",
			low_qual="{folder}/03_filtered/01_bam_filter_low_qual/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
		params:
			q=lambda wildcards: int(get_from_sample_file("MAPQ", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq))
		resources:
			memory=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
			runtime=lambda wildcards, attempt: get_runtime_alloc("filtering", attempt, 24)
		log:
			"{folder}/03_filtered/01_bam_filter/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
		threads: 
			get_threads("filtering", 4)
		conda:
			"../envs/samtools.yaml"
		envmodules:
			module_samtools
		message: "--- SAMTOOLS FILTER {input}"
		shell:
			"""
			samtools view -b --threads {threads} -F 4 -q {params.q} \
			-U {output.low_qual} {input} > {output.mapped} 2> {log}
			"""
else:
	rule samtools_filter:
		"""
		Filter mappings following quality and discard low quality mappings
		"""
		input:
			"{folder}/02_mapped/03_bam_sort/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
		output:
			mapped="{folder}/03_filtered/01_bam_filter/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam",
		params:
			q=lambda wildcards: int(get_from_sample_file("MAPQ", wildcards.id_sample, wildcards.id_library, wildcards.id_fastq))
		resources:
			memory=lambda wildcards, attempt: get_memory_alloc("filtering", attempt, 4),
			runtime=lambda wildcards, attempt: get_runtime_alloc("filtering", attempt, 24)
		log:
			"{folder}/03_filtered/01_bam_filter/{id_sample}/{id_library}/{id_fastq}.{id_genome}.log"
		threads: 
			get_threads("filtering", 4)
		conda:
			"../envs/samtools.yaml"
		envmodules:
			module_samtools
		message: "--- SAMTOOLS FILTER {input}"
		shell:
			"""
			samtools view -b --threads {threads} -F 4 -q {params.q} {input} > {output.mapped} 2> {log}
			"""


##########################################################################################
rule get_final_fastq:
	"""
	Get the final bam file from the fastq part
	"""
	input:
		get_final_bam_fastq
	output:
		"{folder}/04_final_fastq/01_bam/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
	message: "--- GET FINAL BAM {input} (FASTQ LEVEL)"
	run:
		symlink_rev(input, output)
	
rule get_final_fastq_low_qual:
	"""
	Get the final bam file from the fastq part
	"""
	input:
		"{folder}/03_filtered/01_bam_filter_low_qual/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"
	output:
		"{folder}/04_final_fastq/01_bam_low_qual/{id_sample}/{id_library}/{id_fastq}.{id_genome}.bam"	
	message: "--- GET FINAL LOW_QUAL BAM {input} (FASTQ LEVEL)"
	run:
		symlink_rev(input, output)
