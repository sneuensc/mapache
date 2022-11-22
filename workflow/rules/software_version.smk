
## get the version of actually used tools


rule version_snakemake:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/snakemake.txt"
    shell:
        """
        echo snakemake version $(snakemake --version) > {output} || echo snakemake not available > {output}
        """

rule version_adapterremoval:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/adapterremoval.txt"
    conda: "../envs/adapterremoval.yaml"
    envmodules: module_adapterremoval
    shell:
        """
        AdapterRemoval --version 2> >(sed 's/ver./version/g') > {output} || echo AdapterRemoval not available > {output}
        """

rule version_bwa:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/bwa.txt"
    conda: "../envs/bwa.yaml"
    envmodules: module_bwa
    shell:
        """
        set +e;
        txt=$(bwa 2>&1 | grep Version | sed 's/Version:/bwa version/g')
        if [[ "$txt" == *"bwa"* ]]; then
            echo $txt > {output};
        else
            echo bwa not available > {output};
        fi
        exit 0;
        """

rule version_bowtie2:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/bowtie2.txt"
    conda: "../envs/bowtie2.yaml"
    envmodules: module_bowtie2
    shell:
        """
        set +e;
        txt=$(bowtie2 2>&1 | grep "Bowtie 2 version")
        if [[ "$txt" == *"Bowtie"* ]]; then
            echo $txt > {output};
        else
            echo bowtie2 not available > {output};
        fi
        exit 0;
        """

rule version_fastqc:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/fastqc.txt"
    conda: "../envs/fastqc.yaml"
    envmodules: module_fastqc
    shell:
        """
        fastqc --version | sed 's/v/version /g' > {output} || echo fastqc not available > {output}
        """

rule version_gatk3:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/gatk.txt"
    conda: "../envs/gatk3.yaml"
    envmodules: module_gatk3
    shell:
        """
        echo GenomeAnalysisTK version $(GenomeAnalysisTK --version) > {output} || echo GenomeAnalysisTK not available > {output}
        """

rule version_mapdamage2:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/mapdamage2.txt"
    conda: "../envs/mapdamage.yaml"
    envmodules: module_mapdamage
    shell:
        """
        echo mapDamage version $(mapDamage --version) > {output} || echo mapDamage not available > {output}
        """

rule version_picard:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/picard.txt"
    conda: "../envs/picard.yaml"
    envmodules: module_picard
    params: PICARD=get_picard_bin()
    shell:
        """
        set +e;
        txt=$({params.PICARD} MarkDuplicates --version 2>&1 | sed 's/Version:/Picard MarkDuplicates version /g')
        if [[ "$txt" == *"Picard"* ]]; then
            echo $txt > {output};
        else
            echo Picard MarkDuplicates not available > {output};
        fi
        exit 0;
        """

rule version_samtools:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/samtools.txt"
    conda: "../envs/samtools.yaml"
    envmodules: module_samtools
    shell:
        """
        samtools --version | grep samtools | sed 's/samtools/samtools version/g' > {output} || echo samtools not available > {output}
        """

rule version_bcftools:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/bcftools.txt"
    conda: "../envs/bcftools.yaml"
    envmodules: module_bcftools
    shell:
        """
        bcftools --version |grep bcftools | sed 's/bcftools/bcftools version/g' > {output} || echo bcftools not available > {output}
        """

rule version_bamutil:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/bamutil.txt"
    conda: "../envs/bamutil.yaml"
    envmodules: module_bamutil
    shell:
        """
        set +e;
        txt=$(bam 2>&1 | grep Version | sed 's/Version:/BamUtil version/g')
        if [[ "$txt" == *"BamUtil"* ]]; then
            echo $txt > {output};
        else
            echo BamUtil not available > {output};
        fi
        exit 0;
        """

rule version_bedtools:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/bedtools.txt"
    conda: "../envs/bedtools.yaml"
    envmodules: module_bedtools
    shell:
        """
        bedtools --version | sed 's/v/version /g' > {output} || echo bedtools not available > {output}
        """

rule version_seqtk:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/seqtk.txt"
    conda: "../envs/seqtk.yaml"
    envmodules: module_seqtk
    shell:
        """
        set +e;
        txt=$(seqtk 2>&1 | grep Version | sed 's/Version:/seqt version/g')
        if [[ "$txt" == *"seqtk"* ]]; then
            echo $txt > {output};
        else
            echo seqtk not available > {output};
        fi
        exit 0;
        """

rule version_dedup:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/dedup.txt"
    conda: "../envs/dedup.yaml"
    envmodules: module_dedup
    shell:
        """
        dedup --version &> >(head -n1 | sed 's/v/version /g') > {output} || echo DeDup not available > {output}
        """

rule version_qualimap:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/qualimap.txt"
    conda: "../envs/qualimap.yaml"
    envmodules: module_qualimap
    shell:
        """
        qualimap --help | grep QualiMap | sed 's/v./version /g' > {output} || echo qualimap not available > {output}
        """

rule version_multiqc:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/multiqc.txt"
    conda: "../envs/multiqc.yaml"
    envmodules: module_multiqc
    shell:
        """
        multiqc --version | sed 's/,//g' > {output} || echo multiqc not available > {output}
        """

rule version_glimpse:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/glimpse.txt"
    conda: "../envs/glimpse.yaml"
    envmodules: module_glimpse
    shell:
        """
        GLIMPSE_chunk --help | grep Version | sed 's/  \* Version       :/glimpse version/g' > {output} || echo glimpse not available > {output}
        """

rule version_fastp:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/fastp.txt"
    conda: "../envs/fastp.yaml"
    envmodules: module_fastp
    shell:
        """
        fastp --version 2> >(sed 's/fastp/fastp version/g') > {output} || echo fastp not available > {output}
        """

rule version_r:
    output: f"{RESULT_DIR}/04_stats/02_separate_tables/software/R.txt"
    conda: "../envs/r.yaml"
    envmodules: module_r
    shell:
        """
        R --version | head -n 1 > {output} || echo R not available > {output}
        """


