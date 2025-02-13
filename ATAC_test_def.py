import os
import pandas as pd

configfile: "/depot/bioinfo/data/Personal_Directories/arora/atac_test/config.yaml"

shell.prefix(config["module"])

def get_samples():
    input_dir = config["input_dir"]
    samples = sorted(
        set(f.split("_1.fastq.gz")[0] for f in os.listdir(input_dir) if f.endswith("_1.fastq.gz"))
    )
    if not samples:
        raise ValueError("No samples found. Check naming of FASTQ files.")
    print("Detected samples:", samples)
    return samples

SAMPLES = get_samples()

trimmed_dir = config.get("trimmed_dir", "trimmed").rstrip("/")
fastqc_trimmed_dir = config.get("fastqc_trimmed_dir", "fastqc_trimmed").rstrip("/")
bowtie2_index_prefix = config.get("bowtie2_index_prefix")
if not bowtie2_index_prefix:
    raise ValueError("Missing 'bowtie2_index_prefix' in config.")

genome_size = config.get("genome_size")
if not genome_size:
    raise ValueError("Missing 'genome_size' in config.")

for directory in [
    config["fastqc_raw"],
    config["raw_bam"],
    config["sorted_bam"],
    config["clean_bam"],
    config["bedfiles"],
    config["peakcalling"],
    fastqc_trimmed_dir,    
    trimmed_dir
]:
    os.makedirs(directory, exist_ok=True)


###############################################################################
# rule all ##Boiler UP
###############################################################################
rule all:
    input:
        expand(os.path.join(config["bedfiles"], "{sample}.bed"), sample=SAMPLES),
        expand(os.path.join(config["peakcalling"], "{sample}_peaks.narrowPeak"), sample=SAMPLES)

rule index_genome:
    """Generate bowtie2 index files with prefix 'bowtie2_index_prefix'"""
    input:
        config["genome_fa"]
    output:
        # Bowtie2 builds files like prefix.1.bt2, prefix.2.bt2, etc.
        [f"{bowtie2_index_prefix}.{ext}" for ext in 
         ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]]
    shell:
        """
        if [ ! -f {input} ]; then
            echo "Genome FASTA file {input} not found!" >&2
            exit 1
        fi
        bowtie2-build {input} {bowtie2_index_prefix}
        """


rule trim_galore:
    input:
        r1=lambda wc: os.path.join(config["input_dir"], f"{wc.sample}_1.fastq.gz"),
        r2=lambda wc: os.path.join(config["input_dir"], f"{wc.sample}_2.fastq.gz")
    output:
        r1_trimmed=os.path.join(trimmed_dir, "{sample}_1_trimmed.fq.gz"),
        r2_trimmed=os.path.join(trimmed_dir, "{sample}_2_trimmed.fq.gz")
    log:
        os.path.join(trimmed_dir, "{sample}_trim_galore.log")
    shell:
        """
        trim_galore --paired {input.r1} {input.r2} -o {trimmed_dir} > {log} 2>&1
        mv {trimmed_dir}/{wildcards.sample}_1_val_1.fq.gz {output.r1_trimmed}
        mv {trimmed_dir}/{wildcards.sample}_2_val_2.fq.gz {output.r2_trimmed}
        """

rule fastqc_trimmed:
    """Run FastQC on trimmed reads."""
    input:
        r1=os.path.join(trimmed_dir, "{sample}_1_trimmed.fq.gz"),
        r2=os.path.join(trimmed_dir, "{sample}_2_trimmed.fq.gz")
    output:
        html1=os.path.join(fastqc_trimmed_dir, "{sample}_1_fastqc.html"),
        zip1=os.path.join(fastqc_trimmed_dir, "{sample}_1_fastqc.zip"),
        html2=os.path.join(fastqc_trimmed_dir, "{sample}_2_fastqc.html"),
        zip2=os.path.join(fastqc_trimmed_dir, "{sample}_2_fastqc.zip")
    log:
        os.path.join(fastqc_trimmed_dir, "{sample}_fastqc_trimmed.log")
    shell:
        """
        fastqc -o {fastqc_trimmed_dir} {input.r1} {input.r2} > {log} 2>&1
        """


rule align_reads:
    input:
        r1=os.path.join(trimmed_dir, "{sample}_1_trimmed.fq.gz"),
        r2=os.path.join(trimmed_dir, "{sample}_2_trimmed.fq.gz")
    output:
        os.path.join(config["raw_bam"], "{sample}.bam")
    log:
        os.path.join(config["raw_bam"], "{sample}.align_reads.log")
    shell:
        """
        bowtie2 -x {bowtie2_index_prefix} -1 {input.r1} -2 {input.r2} 2> {log} | \
        samtools view -bS - > {output}
        """

rule sort_bam:
    input:
        os.path.join(config["raw_bam"], "{sample}.bam")
    output:
        os.path.join(config["sorted_bam"], "{sample}.sorted.bam")
    log:
        os.path.join(config["sorted_bam"], "{sample}.sort_bam.log")
    shell:
        "samtools sort {input} -o {output} > {log} 2>&1"


rule remove_duplicates:
    input:
        sorted_bam=os.path.join(config["sorted_bam"], "{sample}.sorted.bam")
    output:
        bam=os.path.join(config["clean_bam"], "{sample}.clean.bam"),
        bai=os.path.join(config["clean_bam"], "{sample}.clean.bam.bai")
    log:
        os.path.join(config["clean_bam"], "{sample}.remove_duplicates.log")
    resources:
        picard=1
    shell:
        """
        java -Xmx100g -jar /depot/bioinfo/apps/apps/picard-tools-3.0.0/picard.jar MarkDuplicates \
        I={input.sorted_bam} O={output.bam} M={output.bam}.metrics > {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1
        """

rule convert_to_bed:
    input:
        os.path.join(config["clean_bam"], "{sample}.clean.bam")
    output:
        os.path.join(config["bedfiles"], "{sample}.bed")
    log:
        os.path.join(config["bedfiles"], "{sample}.convert_to_bed.log")
    shell:
        "bedtools bamtobed -i {input} > {output} 2> {log}"


rule all_deduplicated:
    input:
        expand(os.path.join(config["clean_bam"], "{sample}.clean.bam"), sample=SAMPLES)
    output:
        "all_dedup_done.txt"
    shell:
        "touch {output}"


rule call_peaks:
    input:
        bam=os.path.join(config["clean_bam"], "{sample}.clean.bam"),
        dedup_done="all_dedup_done.txt"
    output:
        narrowPeak=os.path.join(config["peakcalling"], "{sample}_peaks.narrowPeak")
    log:
        os.path.join(config["peakcalling"], "{sample}.call_peaks.log")
    shell:
        """
        macs3 callpeak -t {input.bam} \
        -f BAMPE -g {genome_size} -n {wildcards.sample} \
        -q 0.01 --nomodel --shift -75 --extsize 150 --call-summits \
        --keep-dup all --outdir {config[peakcalling]} > {log} 2>&1

        # If you also want broad peaks in the same run:
        macs3 callpeak -t {input.bam} \
        -f BAMPE -g {genome_size} -n {wildcards.sample} \
        -q 0.01 --nomodel --shift -75 --extsize 150 \
        --keep-dup all --broad --outdir {config[peakcalling]} >> {log} 2>&1
        """