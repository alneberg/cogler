from collections import defaultdict

def read_blast_output(blastoutfile): 
    sseq_ids = []
    records = []
    with open(blastoutfile) as in_handle:
        for line in in_handle:
            line_items = line.split("\t")
            qseq = line_items[0]
            sseq = line_items[1]
            pident = line_items[3]
            send = line_items[8]
            sstart = line_items[7]
            slen = line_items[10]

            records.append({'qseqid': qseq,
                            'sseqid': sseq,
                            'pident': float(pident),
                            'send': float(send),
                            'sstart': float(sstart),
                            'slen': float(slen)})

            sseq_ids.append(sseq.split('|')[2])
    return records, sseq_ids

def read_gff_file(gfffile):
    featureid_locations={}
    limits=dict(gff_type=["gene","mRNA","CDS"])
    with open(gfffile) as in_handle:
        for rec in GFF.parse(in_handle, limit_info=limits):
            for feature in rec.features:
                featureid_locations[feature.id] = rec.id
    return featureid_locations

def read_markers_file(marker_file):
    # Stores each line of marker_file as an item in a list
    with open(marker_file) as mf:
        return [l.strip() for l in mf.readlines()]

def read_clustering_file(cluster_file):
    # Returns the cluster names and the contig names per cluster
    contigs_per_cluster = defaultdict(list)
    clusters = set()
    with open(cluster_file) as cf:
        for line in cf.readlines():
            line_items = line.strip().split(',')
            cluster = line_items[1]
            contig = line_items[0]
            clusters.add(cluster)
            contigs_per_cluster[cluster].append(contig)
    return list(clusters), contigs_per_cluster


