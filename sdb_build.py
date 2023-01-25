import os

otus = "/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/vog_transcript_otus/"
outfile = "/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/vog_530_3.tsv"
viral_tax = "/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/viral_taxonomy_fixed.tsv"

id_to_tax = {}
with open(viral_tax) as r:
    for line in r:
        arr = line.split("\t")
        id = arr[0]
        tax = "; ".join(arr[1][:-1].split(";"))
        id_to_tax[arr[0]] = tax

id_to_tax

with open(outfile, "w+") as w:
    w.write("gene\tsample\tsequence\tnum_hits\tcoverage\ttaxonomy\tread_names\tnucleotides_aligned\ttaxonomy_by_known?\tread_unaligned_sequences\tequal_best_hit_taxonomies\ttaxonomy_assignment_method\n")
    for otu in os.listdir(otus):
        count = 0
        with open(otus+otu) as r:
            for line in r.readlines()[1:]:
                arr = line.split('\t')
                if arr[0] == arr[1]:
                    if len(arr[6].split(" ")) > 1:
                        continue
                    arr[5] = id_to_tax[arr[6]]
                    w.write(line)
                    count += 1
        if count == 0:
            print("No correctly-aligned OTUs found for: {}".format(otu))
