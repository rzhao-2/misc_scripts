import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    matches_fasta, mfqe, output, hmms_and_names, taxfiles,  log = params
    logging.debug("Processing genome: " + genome)

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = f"python scripts/transpose_hmms_with_sequences.py --input-fasta {matches_fasta} --taxonomy {' '.join(taxfiles)} --hmm-seq {mfqe} --hmm-spkg {hmms_and_names} --output {output} &> {log}"
    # logging.info(cmd)
    extern.run(cmd)

protein_filepaths = [filepath.strip() for filepath in open(snakemake.params.protein_filepaths)]
matches_dir = snakemake.params.matches_dir
mfqe_dir = snakemake.params.mfqe_dir
taxfiles = snakemake.params.taxfiles
output_dir = snakemake.params.output_dir

hmms_and_names = snakemake.params.hmms_and_names

logs_dir = snakemake.params.logs_dir
num_threads = snakemake.threads

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(protein_filepaths), num_threads))

pathlib.Path(os.path.dirname(output_dir)).mkdir(parents=True, exist_ok=True)

param_sets = []
for filepath in protein_filepaths:
    genome = os.path.basename(filepath).rsplit(".", 1)[0]
    # fasta = output_dir + "/hmmsearch/matches/{genome}",
    # matches = output_dir + "/hmmsearch/matches/{genome}.fam",
    matches_fasta = os.path.join(mfqe_dir, genome + '.faa')
    fam = os.path.join(matches_dir, genome + ".tsv")
    output = os.path.join(output_dir, genome)
    log = os.path.join(logs_dir, f"{genome}_transpose.log")

    if not os.path.exists(fam):
        logging.warning(f"Missing fam file for {genome}")
        continue
    # pfam_search, tigrfam_search, hmms_and_names, output, log
    param_sets.append((matches_fasta, fam, output, hmms_and_names, taxfiles, log))

process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass