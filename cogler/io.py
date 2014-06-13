from collections import defaultdict
import pandas as pd

def read_per_phylum_scgs(fn):
    """Read information from phylum scg file.

    Arguments:
    fn -- file name

    Output:
    summary -- pandas dataframe with information for each phylum
    phyla_scgs -- pandas dataframe with cogs for each phylum
    """
    df = pd.read_table(fn, index_col=0)
    summary = df.ix[:,'Number_genomes':'Number_SCG']
    phyla_scg = df.ix[:, 'COG0001':]
    return summary, phyla_scg

def read_blast_output(blastoutfile): 
    """Parses rpsblast output file.

    Arguments:
    blastoutfile -- file name

    Output:
    records -- list of rpsblast output lines, each line is a dictionary
               with the column headers as keys and the column values
               as values.
    sseq_ids -- list of subject ids, cdd ids in this case.

    Assumes rpsblast has been run with the -outfmt "6 qseqid sseqid
    evalue pident score qstart qend sstart send length slen" option.
    """
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
    """Read file with list of marker gene names.

    Arguments:
    marker_file -- file name

    Output:
        -- List of all marker genes.
    """
    with open(marker_file) as mf:
        return [l.strip() for l in mf.readlines()]

def read_clustering_file(cluster_file):
    """Read CONCOCT style clustering file, with contig_id, cluster_id pairs.

    Arguments:
    cluster_file -- file name

    Output:
    cluster_per_contig -- dictionary with contig ids as keys and cluster ids as
                           values.
    """
    with open(cluster_file, 'r') as cf:
        cluster_per_contig = dict([tuple(row.strip().split(',')) for row in cf.readlines()])
    return cluster_per_contig

