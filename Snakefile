import glob

configfile: "config.yaml"


rule all:
    input: f"{config['output_dir']}/presence_absence_matrix.txt"


if config["BAKTA"]["exec"] and config["GFF"]["exec"]:
    print("cannot have BAKTA and GFF true in config file")
    exit(1)


# if wanting to annotate MAGs with BAKTA
if config["BAKTA"]["exec"]:
    rule bakta:
        input:
            genome = f"{config['genome_fasta']}/{{sample}}.fasta"
        output:
            ann_dir = directory(f"{config['output_dir']}/annotated/{{sample}}_ann")
        conda:
            "bakta"
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 15000
        params:
            trn = config['training_file'],
            DB = directory(f"{config['bakta_db']}")
        log:
            f"{config['output_dir']}/logs/bakta/{{sample}}.log"
        shell:
            """
            bakta {input.genome} --db {params.DB} --prefix {wildcards.sample} --prodigal-tf {params.trn} \
             --translation-table 11 --threads {threads} --output {output.ann_dir} >{log} 2>&1
            """


    def get_samples(genome_dir):
        list_of_samples = glob.glob(genome_dir + "/*.fasta")
        new_list_of_samples = []
        for sam in list_of_samples:
            sam = sam.split("/")[-1]
            sam = sam.replace(".fasta", "")
            new_list_of_samples.append(sam)
        return new_list_of_samples


    rule fix_ffn_file:
        input:
            annotations = expand(f"{config['output_dir']}/annotated/{{sample}}_ann", sample=get_samples(config['genome_fasta']))
        output:
            fixed_annotations = directory(f"{config['output_dir']}/all_ffn")
        conda: #just needs biopython
            "Snakemake"
        script: "scripts/fix_ffn_files.py"


# if you have existing GFF to use:
elif config["GFF"]["exec"]:
    rule transform_gff:
        input:
            gff_dir = f"{config['gff_dir']}"
        output:
            ffn_files = directory(f"{config['output_dir']}/all_ffn")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 15000
        shell:
            """
            python scripts/panphlan_pangenome_generation.py --i_gff {input}
            mv ffn_from_gff {output.ffn_files}
            """


rule concat:
    input:
        ffn_files = f"{config['output_dir']}/all_ffn"
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
        f"{config['output_dir']}/logs/mmseqs2.log"
    params:
        seq_id = 0.9,
        cov_mode = 0,
        c = 0.8,
        output_prefix = f"{config['output_dir']}/mmseqs/mmseqs",
        tmp_dir = f"{config['output_dir']}/mmseqs/tmp"
    conda:
        "mmseq2"
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


rule sort_mmseqs2_clusters:
    input:
        clusters = rules.mmseqs2.output.clusters
    output:
        sorted_clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.sorted.tsv"
    threads: 1
    resources:
        mem_mb=2000
    log:
        f"{config['output_dir']}/logs/sort_mmseqs2_clusters.log"
    shell: "sort {input.clusters} > {output.sorted_clusters}"


rule build_matrix:
    input:
        rep_list = f"{config['output_dir']}/mmseqs/rep_sequences.list",
        clusters = f"{config['output_dir']}/mmseqs/mmseqs_cluster.sorted.tsv"
    output:
        matrix = f"{config['output_dir']}/presence_absence_matrix.txt"
    threads: 1
    resources:
        mem_mb=5000
    conda:
        # env has biopython and pandas
        "poppunk"
    script: "scripts/make_presence_absence_matrix.py"

