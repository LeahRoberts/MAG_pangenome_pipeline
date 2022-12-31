import glob
from Bio import SeqIO


# inputdir and outputdir need to be without trailing "/"

def fix_ffn_files(annotations_dir, output_dir):
    fastas = glob.glob(annotations_dir + "/*/*ffn")
    for file in fastas:
        filename = file.split("/")[-1]
        filename = filename.replace(".ffn", "")
        output_handle = output_dir + "/" + filename + ".formatted.ffn"
        records = SeqIO.parse(file, "fasta")
        all_records = []
        for record in records:
            short_id = record.id.split("_")[1]
            record.id = filename + "_" + short_id
            all_records.append(record)
        SeqIO.write(all_records, output_handle, "fasta")


fix_ffn_files(snakemake.input.annotations, snakemake.output.fixed_annotations)
