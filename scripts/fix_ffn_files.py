import glob
from Bio import SeqIO
import os
# inputdir and outputdir need to be without trailing "/"


def fix_ffn_files(list_of_annotation_dirs, outputdir):
    try:
        os.mkdir(outputdir)
    except OSError:
        pass
    for ann_dirs in list_of_annotation_dirs:
        file = glob.glob(ann_dirs + "/*.ffn")[0] 
        filename = file.split("/")[-1]
        filename = filename.replace(".ffn", "")
        records = SeqIO.parse(file, "fasta")
        outfilename = outputdir + "/" + filename + "_fixed.ffn"
        all_records = []
        for record in records:
            short_id = record.id.split("_")[1]
            record.id = filename + "_" + short_id
            all_records.append(record)
        SeqIO.write(all_records, outfilename, "fasta")


fix_ffn_files(snakemake.input.annotations, snakemake.output.fixed_annotations)
