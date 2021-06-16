
##########################################################################################
## all rules for the stats
##########################################################################################
## individual stats

rule fastqc:
    """
    Quality control of fastq file by fastqc (SE or R1)
    """
    input:
        "results/{folder}/{file}.fastq.gz"
    output:
        html="results/{folder}/{file}_fastqc.html",
        zip="results/{folder}/{file}_fastqc.zip"
    log:
        "results/logs/{folder}/{file}_fastqc.log"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("fastqc_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("fastqc_time", attempt, 1)
    conda:
    	"../envs/fastqc.yaml"
    envmodules:
    	module_fastqc
    message: "--- FASTQC {input}"
    shell:
        "fastqc --quiet --outdir $(dirname {output.html}) {input} 2> {log}"


rule samtools_flagstat:
    """
    Compute samtools falgstat on bam file
    """
    input:
        bam="results/{folder}/{file}.bam"
    output:
        "results/{folder}/{file}_flagstat.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("samtools_flagstat_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("samtools_flagstat_time", attempt, 1)
    params: ""
    log:
        "results/logs/{folder}/{file}_flagstat.log"
    conda:
        "../envs/samtools.yaml"
    envmodules:
        module_samtools
    message: "--- SAMTOOLS FLAGSTAT {input}"        
    shell:
        "samtools flagstat --threads {threads} {input} > {output} 2> {log};"

##########################################################################################
## merging individual stats
rule multiqc_fastqc:
    """
    Merging fastqc output with multiqc
    """
    input:
        fastq = get_fastqc_for_multiqc
    output:
        #html = report("results/{folder}/multiqc_fastqc_{group}.html", category=" Quality control"),
        html = "results/{folder}/multiqc_fastqc_{group}.html",
        txt = "results/{folder}/multiqc_fastqc_{group}_data/multiqc_fastqc.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("multiqc_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("multiqc_time", attempt, 1)
    log:
        "results/logs/{folder}/multiqc_fastqc_{group}.log"
    conda:
        "../envs/multiqc.yaml"
    message: "--- MULTIQC fastqc_{wildcards.group}: {input}"
    shell:
        """
        multiqc -n $(basename {output.html}) -f -d -o $(dirname {output.html}) {input} \
        -b \"Merged fastqc output of {wildcards.group} fastq files\" 2> {log}
        """
        

rule multiqc_flagstat:
    """
    Merging fastqc output with multiqc
    """
    input:
        get_flagstat_for_multiqc
    output:
        html = "results/{folder1}/{folder2}/{folder3}/multiqc_flagstat_{group}.{id_genome}.html",
        txt = "results/{folder1}/{folder2}/{folder3}/multiqc_flagstat_{group}.{id_genome}_data/multiqc_samtools_flagstat.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("multiqc_mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("multiqc_time", attempt, 1)
    log:
        "results/logs/{folder1}/{folder2}/{folder3}/multiqc_flagstat_{group}.{id_genome}.log"
    conda:
        "../envs/multiqc.yaml"
    message: "--- MULTIQC flagstat_{wildcards.group}: {input}"
    shell:
        """
        multiqc -n $(basename {output.html}) -f -d -o $(dirname {output.html}) {input} \
        -b \"Merged flagstat output of {wildcards.group} libraries and genome {wildcards.id_genome}\" 2> {log}
        """

## write summary stats for fastq
rule write_summary_stats_fastq:
    """
    Combine fastq statistics across files and write them to file
    """
    input:
        orig = "results/01_fastq/05_stats/01_multiqc/multiqc_fastqc_orig_data/multiqc_fastqc.txt",
        trim = "results/01_fastq/05_stats/01_multiqc/multiqc_fastqc_trim_data/multiqc_fastqc.txt",
        final = "results/01_fastq/05_stats/01_multiqc/multiqc_flagstat_final.{id_genome}_data/multiqc_samtools_flagstat.txt"
    output: 
        fastq_stats = report("results/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.csv", caption="../report/fastq_stats.rst", category="Mapping statistics table")
    log: 
        "results/logs/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.log"
    message: "--- WRITE FASTQ SUMMARY STATISTICS OF {wildcards.id_genome}"
    run:
        import pandas as pd
        import numpy as np
        
        def read_flagstat(file, extract, name):
            d1 = pd.read_csv(file, sep="\t").rename({extract: name}, axis=1)[['Sample', name]] 
            d1[name] = d1[name].astype('int64')
            d2 = pd.DataFrame(d1['Sample'].str.split('|').tolist(), columns = ['type','type','type','type','SM','LB','ID'])
            d2['SM'] = d2['SM'].str.strip()
            d2['LB'] = d2['LB'].str.strip()
            d2['ID'] = d2['ID'].str.strip()   
            if "_R1_data" == os.path.dirname(file)[-8:]:
                d2['ID'] = d2['ID'].map(lambda x: str(x)[:-3])
            d2['ID'] = d2['ID'].map(lambda x: x.rstrip(".{}_flagstat".format(wildcards.id_genome)))
            dd = pd.concat([d2[['SM','LB','ID']], d1[name]], axis=1)
            return (dd)

        ## sample file    
        db_fastq = pd.read_csv(SAMPLES, sep=delim)[['SM', 'LB' ,'ID']]\
        .sort_values(['SM', 'LB' ,'ID'], axis=0, ascending=True)\
        .reset_index(drop=True)
        
        ## fastq file
        db_orig = read_flagstat(input.orig, 'Total Sequences', 'reads_raw')
        db_trim = read_flagstat(input.trim, 'Total Sequences', 'reads_trim')
        db_map = read_flagstat(input.final, 'mapped_passed', 'mapping')
        
        db_fastq = pd.merge(db_fastq, db_orig, how="left", on = ['SM','LB','ID'])
        db_fastq = pd.merge(db_fastq, db_trim, how="left", on = ['SM','LB','ID'])
        db_fastq = pd.merge(db_fastq, db_map, how="left", on = ['SM','LB','ID'])
        db_fastq['trim_prop'] = np.where(db_fastq['reads_raw'].ne(0), db_fastq['reads_trim'].astype(float) / db_fastq['reads_raw'], np.NaN) 
        db_fastq['endo_prop'] = np.where(db_fastq['reads_raw'].ne(0), db_fastq['mapping'].astype(float) / db_fastq['reads_raw'], np.NaN) 
        
        db_fastq[['ID', 'LB', 'SM', 'reads_raw', 'reads_trim', 'trim_prop', 'mapping', 'endo_prop']]\
        .round({'reads_trim': 0, 'trim_prop': 4, 'mapping': 0, 'endo_prop': 4})\
        .to_csv(output[0], index = None, header=True)


## write summary stats for libraries
rule write_summary_stats_library:
    """
    Combine library statistics across files and write them to file
    """
    input:
        fastq = "results/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.csv",
        final = "results/02_library/04_stats/01_multiqc/multiqc_flagstat_final.{id_genome}_data/multiqc_samtools_flagstat.txt"
    output: 
       fastq_stats = report("results/02_library/04_stats/02_summary/library_stats.{id_genome}.csv", caption="../report/library_stats.rst", category="Mapping statistics table")
    log: 
        "results/logs/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.log"
    message: "--- WRITE LIBRARY SUMMARY STATISTICS OF {wildcards.id_genome}"
    run:
        import pandas as pd
        import numpy as np
        
        def read_flagstat(file, extract, name):
            d1 = pd.read_csv(file, sep="\t").rename({extract: name}, axis=1)[['Sample', name]] 
            d1[name] = d1[name].astype('int64')
            d2 = pd.DataFrame(d1['Sample'].str.split('|').tolist(), columns = ['type','type','type','type','SM','LB'])
            d2['SM'] = d2['SM'].str.strip()
            d2['LB'] = d2['LB'].str.strip()
            d2['LB'] = d2['LB'].map(lambda x: x.rstrip(".{}_flagstat".format(wildcards.id_genome)))
            dd = pd.concat([d2[['SM','LB']], d1[name]], axis=1)
            return (dd)

        ## fastq stats
        db_fastq = pd.read_csv(input.fastq)
        db_library = db_fastq.groupby(['SM','LB'])[['reads_raw', 'reads_trim', 'mapping']].apply(lambda x: x.sum(skipna=False)).reset_index()
        
        ## library stats
        db_lib = read_flagstat(input.final, 'mapped_passed', 'mapping_final')
        
        db_library = pd.merge(db_library, db_lib, how="left", on = ['SM','LB'])
        db_library['duplicates'] = db_library['mapping'] - db_library['mapping_final']
        db_library['trim_prop'] = np.where(db_library['reads_raw'].ne(0), db_library['reads_trim'].astype(float) / db_library['reads_raw'], np.NaN) 
        db_library['endo_prop'] = np.where(db_library['reads_raw'].ne(0), db_library['mapping'].astype(float) / db_library['reads_raw'], np.NaN) 
        db_library['endo_final_prop'] = np.where(db_library['reads_raw'].ne(0), db_library['mapping_final'].astype(float) / db_library['reads_raw'], np.NaN) 
        db_library['duplicates_prop'] = np.where(db_library['mapping'].ne(0), db_library["duplicates"].astype(float) / db_library["mapping"], np.NaN) 

        db_library[['LB', 'SM','reads_raw', 'reads_trim', 'trim_prop', 'mapping', 'endo_prop', 'duplicates', 'duplicates_prop', 'mapping_final', 'endo_final_prop']]\
        .round({'reads_raw': 0, 'reads_trim': 0, 'trim_prop': 4, 'mapping': 0, 'endo_prop': 4, 'duplicates': 0, 'duplicates_prop': 4, 'mapping_final': 0, 'endo_final_prop': 4})\
        .to_csv(output[0], index = None, header=True)


## write summary stats for sample
rule write_summary_stats_sample:
    """
    Combine sample statistics across files and write them to file
    """
    input:
        library = "results/02_library/04_stats/02_summary/library_stats.{id_genome}.csv"
    output: 
        sample_stats = report("results/03_sample/04_stats/01_summary/sample_stats.{id_genome}.csv", caption="../report/sample_stats.rst", category="Mapping statistics table")
    log: 
        "results/logs/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.log"
    message: "--- WRITE SAMPLE SUMMARY STATISTICS OF {wildcards.id_genome}"
    run:
        import pandas as pd
        import numpy as np
        
        ## sample stats
        db_library = pd.read_csv(input.library)
        db_sample = db_library.groupby('SM')[['reads_raw', 'reads_trim', 'mapping', 'mapping_final', 'duplicates']].apply(lambda x: x.sum(skipna=False)).reset_index()
        db_sample['trim_prop'] = np.where(db_sample['reads_raw'].ne(0), db_sample['reads_trim'].astype(float) / db_sample['reads_raw'], np.NaN) 
        db_sample['endo_prop'] = np.where(db_sample['reads_raw'].ne(0), db_sample['mapping'].astype(float) / db_sample['reads_raw'], np.NaN) 
        db_sample['endo_final_prop'] = np.where(db_sample['reads_raw'].ne(0), db_sample['mapping_final'].astype(float) / db_sample['reads_raw'], np.NaN) 
        db_sample['duplicates_prop'] = np.where(db_sample['mapping'].ne(0), db_sample["duplicates"].astype(float) / db_sample["mapping"], np.NaN) 

        db_sample[['SM','reads_raw', 'reads_trim', 'trim_prop', 'mapping', 'endo_prop', 'duplicates', 'duplicates_prop', 'mapping_final', 'endo_final_prop']]\
        .round({'reads_raw': 0, 'reads_trim': 0, 'trim_prop': 4, 'mapping': 0, 'endo_prop': 4, 'duplicates': 0, 'duplicates_prop': 4, 'mapping_final': 0, 'endo_final_prop': 4})\
        .to_csv(output[0], index = None, header=True)


rule plot_summary_statistics:
    """
    Plot summary statistics
    """
    input:
        fastq_stats = "results/01_fastq/05_stats/02_summary/fastq_stats.{id_genome}.csv",
        library_stats = "results/02_library/04_stats/02_summary/library_stats.{id_genome}.csv",
        sample_stats = "results/03_sample/04_stats/01_summary/sample_stats.{id_genome}.csv"
    output: 
        plot_1_nb_reads = report("results/03_sample/04_stats/01_summary/1_nb_reads.{id_genome}.png", caption="../report/1_nb_reads.rst", category="Mapping statistics plots"),
        plot_2_mapped = report("results/03_sample/04_stats/01_summary/2_mapped.{id_genome}.png", caption="../report/2_mapped.rst", category="Mapping statistics plots"),        
        plot_3_endogenous = report("results/03_sample/04_stats/01_summary/3_endogenous.{id_genome}.png", caption="../report/3_endogenous.rst", category="Mapping statistics plots"),        
        plot_4_duplication = report("results/03_sample/04_stats/01_summary/4_duplication.{id_genome}.png", caption="../report/4_duplication.rst", category="Mapping statistics plots")      
    log: 
        "results/logs/03_sample/04_stats/01_summary/plot_summary_statistics_{id_genome}.log"
    conda:
    	"../envs/r.yaml"
    envmodules:
    	module_r
    message: "--- PLOT SUMMARY STATISTICS OF {wildcards.id_genome}"
    script:
    	"../scripts/plot_stats.R"


##########################################################################################
## read depth and sex de3termination
rule depth_compute:
    """
    Compute the read depth per chromosome 
    """
    input:
        "results/{folder}/{file}.{id_genome}.bam"
    output:
        "results/{folder}/{file}.{id_genome}_depth.txt"
    resources:
        memory=lambda wildcards, attempt: get_memory_alloc("mem", attempt, 2),
        runtime=lambda wildcards, attempt: get_runtime_alloc("time", attempt, 1)
    log:
        "results/logs/{folder}/{file}.{id_genome}.log"
    threads: 1
    message: "--- COMPUTE DEPTH OF {input}"
    shell:
        "workflow/scripts/depth.py {input} > {output} 2> {log}"
        
        
rule concat_depth_statistics:
    input:
        get_depth_files
    output:
        report("results/{folder}/depth_stats_{id_genome}.csv", caption="../report/sample_depth.rst", category="Mapping statistics table")
    threads: 1
    log:
        "results/logs/{folder}/depth_stats_{id_genome}.log"
    message: "--- COMBINE DETH STATISTICS TO {output} ---"
    script:
    	"../scripts/depth_concat.py"


rule plot_depth_statistics:
    input:
        sample_depth = "results/{folder}/depth_stats_{id_genome}.csv"
    output:
        plot_5_AvgReadDepth = report("results/{folder}/5_AvgReadDepth.{id_genome}.svg", caption="../report/5_AvgReadDepth.rst", category="Mapping statistics plots"),
        plot_6_AvgReadDepth_MT = report("results/{folder}/6_AvgReadDepth_MT.{id_genome}.svg", caption="../report/6_AvgReadDepth_MT.rst", category="Mapping statistics plots"),
        plot_7_Sex = report("results/{folder}/7_Sex.{id_genome}.svg", caption="../report/7_Sex.rst", category="Mapping statistics plots")
    threads: 1
    log:
        "results/logs/{folder}/depth_stats_plot.{id_genome}.csv.log"
    conda:
    	"../envs/r.yaml"
    envmodules:
    	module_r
    message: "--- PLOT DEPTH STATISTICS OF {input}"
    script:    	
    	"../scripts/plot_depth.R"