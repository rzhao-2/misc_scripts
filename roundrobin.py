import readline

m_thresh = input("input total marker threshold (no input defaults to 3 per species): ")

if m_thresh == '':
    m_thresh = 0
else: 
    m_thresh = int(m_thresh)

def sort_markers(markers):
    sorter = []
    for marker in markers:
        sorter.append((len(sc_markers[marker]), marker))
    sorter.sort(reverse=True)
    sorted_markers = [tup[1] for tup in sorter]
    return sorted_markers

# round robin code
genomes_list = "/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/species-hmmhits-sc.tsv"
marker_list = "/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/vog-hits2.tsv"

sc_markers = {}
species = {}
final_markers = set()
coverage = {}
genome_loop = []

#populate marker dictionary
with open(marker_list) as r:
    for line in r.readlines()[1:]:
        lsplit = line.strip('\n').split('\t')
        vog = lsplit[0]
        genomes = lsplit[5].split(';')
        sc_markers[vog] = genomes

#populate species dictionary
with open(genomes_list) as r:
    for line in r.readlines()[1:]:
        genome, __, __, markers = line.strip('\n').split('\t')
        markers = markers.split(';')
        if len(markers) > 0:
            if markers[0] == "":
                continue
            species[genome] = sort_markers(markers)
            coverage[genome] = 0
            genome_loop.append(genome)

i = 0
while len(final_markers) < m_thresh:
    print("Roundrobin iteration:", i+1)
    for genome in genome_loop:
        if coverage[genome] > i:
            continue
        for marker in species[genome]:
            if marker not in final_markers:
                final_markers.add(marker)
                for covered_genome in sc_markers[marker]:
                    if covered_genome in coverage:
                        coverage[covered_genome] += 1
                break
    i += 1

print(len(final_markers))
with open("/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/picked_markers.txt", "w+") as w:
    w.write('\n'.join([marker for marker in final_markers]))
with open("/mnt/hpccs01/work/microbiome/msingle/rossenzhao/vogdb211/coverage_values.txt", "w+") as w:
    i = 0
    for entry in coverage:
        if coverage[entry] < 15:
            i+=1
        w.write(str(entry)+'\t'+str(coverage[entry])+'\n')
print(i)