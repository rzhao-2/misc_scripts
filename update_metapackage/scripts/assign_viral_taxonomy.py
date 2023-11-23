import os
import sys
import argparse
import logging

import pandas as pd

input_dir = snakemake.params.input_dir
input_files = [f for f in os.listdir(input_dir)]
id_to_tax = {}
with open(snakemake.params.viral_taxonomy) as r:
    for line in r:
        id_, tax = line.strip().split("\t")
        id_to_tax[id_] = tax
#for each input file, read in first line as header, change taxonomy column in all rows to viral taxonomy based on filename, concatenate to single output file
output = pd.DataFrame()
for f in input_files:
    taxonomy = id_to_tax[f.rsplit(".", 2)[0]]
    logging.info(f"Reading {f}")
    df = pd.read_csv(os.path.join(input_dir, f), sep="\t", header=0)
    for row in df.itertuples():
        df.at[row.Index, "taxonomy"] = taxonomy
    output = pd.concat([output, df])

logging.info(f"Writing output to {snakemake.output[0]}")
output.to_csv(snakemake.output[0], sep="\t", index=False)