#1. how to batch submit (for Bakta)
#2. how to wait until one rule finished before starting another
#3. parameter/config to say if already have annotations? i.e. start at rule fix_ffn_file?


configfile: "config.yaml"


rule all:
    input: f"{config['output_dir']}/presence_absence_matrix.txt"

rule bakta:
    input:
        genome = f"{config['genome_fasta']}/{sample}.fasta"
    output:
        ann_dir = directory(f"{config['output_dir']}/annotated"),
        touch(f"{config['output_dir']}/bakta.done")
    conda:
        "envs/bakta.yaml"
    threads: 16
    params:
        trn = config['training_file'],
        output_prefix = "{sample}",
        DB = directory(f"{config['bakta_DB']}")
    log:
        "logs/bakta/{sample}.log"
    shell:
        """
        bakta {input.genome} --db {params.DB} --prefix {params.output_prefix} --prodigal-tf {params.trn} \
         --translation-table 11 --threads {threads} --output {output.ann_dir} >{log} 2>&1
         """


rule fix_ffn_file:
    input:
        "bakta.done",
        annotations = directory(f"{config['output_dir']}/annotated")
    output:
        fixed_annotations = directory(f"{config['output_dir']}/all_ffn")
    conda:
        "envs/biopython.yaml"
    script: "scripts/fix_ffn_files.py"


rule concat:
    input:
        ffn_files = directory(f"{config['output_dir']}/all_ffn")
    output:
        f"{config['output_dir']}/all_samples.concat.ffn"
    shell:
        "cat {input.ffn_files}/* > {output}"


rule mmseqs2:
    input:
        f"{config['output_dir']}/all_samples.concat.ffn"
    output:
        all_seqs = f"{config['output_dir']}/mmseqs/mmseqs_all_seqs.fasta",
        clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.tsv",
        rep_seq = f"{config['output_dir']}/mmseqs/mmseqs_rep_seq.fasta"
    threads: 16
    resources:
        mem_mb=lambda wildcards, attempt: 16000 * attempt
    log:
        "logs/mmseqs2.log"
    params:
        seq_id = 0.9,
        cov_mode = 0,
        c = 0.8,
        output_prefix = f"{config['output_dir']}/mmseqs/mmseqs",
        tmp_dir = f"{config['output_dir']}/mmseqs/tmp"
    conda:
        "envs/mmseqs2.yaml"
    shell:
        "mmseqs easy-cluster {input} {params.output_prefix} {params.tmp_dir} --min-seq-id {params.seq_id} \
        --cov-mode {params.cov_mode} -c {params.c} --threads {threads} >{log} 2>&1"


rule rep_seq_list:
    input:
        f"{config['output_dir']}/mmseqs/mmseqs_rep_seq.fasta"
    output:
        f"{config['output_dir']}/mmseqs/rep_sequences.list"
    shell:
        "grep '>' {input} | cut -f2 -d'>' | cut -f1 -d' ' > {output}"


rule build_matrix:
    input:
        rep_list = f"{config}['output_dir']}/mmseqs/rep_sequences.list",
        clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.tsv"
    output:
        f"{config['output_dir']}/presence_absence_matrix.txt"
    conda:
        # also pandas
        "envs/biopython.yaml"
    script: "scripts/make_presence_absence_matrix.py"

