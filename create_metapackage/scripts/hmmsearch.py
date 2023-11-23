import os

import logging
import pathlib
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    genome_filepath, hmmsearch_output = params  
    logging.debug(f"Processing {genome_filepath}")
    cmd = f"hmmsearch -E 0.00001 --cpu 1 --tblout {hmmsearch_output} {hmm} {genome_filepath} > /dev/null"
    os.system(cmd)

hmm = snakemake.input.hmm
genome_proteins = [prot_filepath.strip('\n') for prot_filepath in open(snakemake.input.genome_proteins)]
# genome_proteins = snakemake.input.protein_list
output_dir = snakemake.params.output_dir
num_threads = snakemake.threads

pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genome protein files with {} threads".format(len(genome_proteins), num_threads))

param_sets = []
for genome_filepath in genome_proteins:
    hmmsearch_output = os.path.join(output_dir, os.path.splitext(os.path.basename(genome_filepath))[0] + ".txt")
    param_sets.append((genome_filepath, hmmsearch_output))
logging.info(os.path.basename(__file__)+ f": Processing {len(param_sets)} genomes with {num_threads} threads")
process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

with open(snakemake.output[0], 'w') as _: pass