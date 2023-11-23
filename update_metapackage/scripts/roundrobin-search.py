import os
import heapq
import logging

# take list of hmmsearch tblout files, iterate through list of genomes and select hmm that covers the most genomes

def get_single_copy_hmms(hmms):
    ratio = 1.05
    kept_hmms = {}
    genome_heapqs = {}
    for hmm in hmms:
        genomes = set()
        num_hits = 0
        # TODO: read hmm tblout, count number of hits per genome
        # if genome has only one hit, add to genomes
        copynumber = num_hits / len(genomes)
        # if number of hits per genome <= ratio 
        if copynumber <= ratio:
            # add hmm to number of kept hmms
            kept_hmms[hmm] = genomes
            for genome in genomes:
                if genome not in genome_heapqs:
                    genome_heapqs[genome] = heapq.heapify([])
                genome_heapqs[genome].heappush((-len(genomes), hmm))
    return kept_hmms, genome_heapqs


def roundrobin(kept_hmms, genome_heapqs):
    picked_hmms = set()
    genomes = {}
    i = 0
    # pick at least 500 HMMs
    while len(picked_hmms) < 500:
        for genome in genome_heapqs:
            if genome not in genomes:
                genomes[genome] = 0
            if genomes[genome] > i:
                continue
            while len(genome_heapqs[genome]) > 0:
                if genome_heapqs[genome][0][1] not in picked_hmms:
                    hmm = genome_heapqs[genome].pop()[1]
                    picked_hmms.add(hmm)
                    for g in kept_hmms[hmm]:
                        if g not in genomes:
                            genomes[g] = 0
                        genomes[g] += 1
                    break
        i+=1
    return picked_hmms

hmms = snakemake.params.hmms
num_packages = snakemake.params.num_packages
output_dir = snakemake.params.output_dir
output_hmm_list = os.path.join(output_dir, "hmm_list.txt")

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} hmms".format(len(hmms)))

kept_hmms, heapqs = get_single_copy_hmms(hmms)
logging.info(os.path.basename(__file__) + ": Found {} single copy hmms for greedy search".format(len(kept_hmms)))

picked_hmms = roundrobin(kept_hmms, heapqs)
logging.info(os.path.basename(__file__) + ": Picked a total of {} hmms from greedy search".format(len(picked_hmms)))
with open(output_hmm_list, "w+") as w:
    for hmm in picked_hmms:
        w.write(hmm + "\n")

logging.info('done')
with open(snakemake.output[0], 'w') as _: pass