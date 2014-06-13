"""Cogler is the COG annotation counter. """
from collections import Counter
import pandas as pd

class CogsPerContig(object):
    """
    Count cogs per contig.

    Constructor Arguments

    ``cluster_per_contig``

        Dictionary with contig ids as keys and cluster id as values.

    ``features_per_contig``

        Dictionary with contig ids as keys and a list of COG accession
        identifiers as values.

    On initialization, a pandas dataframe self.df is created containing
    one row per contig, and one column per cog accession. The values in
    the dataframe is either 1 or NaN, corresponding to if the cog
    was observed in a contig or not.
    """

    def __init__(self, cluster_per_contig, features_per_contig):
        self.cluster_per_contig = cluster_per_contig
        self.features_per_contig = features_per_contig
        self.df = self._construct_count_df()

    def _construct_count_df(self):
        rows = {}
        for key, value in self.features_per_contig.iteritems():
            row = dict(Counter(value))
            row['cluster'] = self.cluster_per_contig[key]
            rows[key] = row
        return pd.DataFrame.from_dict(rows, orient='index')

class Phylum(object):
    """
    Class keeping values corresponding to a phylum.

    Constructor Arguments

    ``name``

        The name of the phylum.

    ``scgs``

        A list of COG accession identifiers found to be single copy genes
        (SCG):s for this phylum.
    """
    def __init__(self, name, scgs):
        self.name = name
        self.scgs = scgs

class Output(object):
    """
    Class to construct the final cog per phyla dataframe.

    Constructor Arguments

    ``cogs_per_contig``

        A CogsPerContig instance.
    """
    def __init__(self, cogs_per_contig):
        self.cogs_per_contig = cogs_per_contig
        self.phyla = []
        self.result_matrix = None

    def add_phylum(self, phylum):
        """
        Add a new phylum to the result matrix, calculates the
        total scgs present in that phylum from the list of scgs in
        the phylum, the number of COGS present in more than 1 copy
        for each cluster and the number of COGS present in exactly
        1 copy.
        """
        assert phylum not in self.phyla
        scgs = filter(lambda c: c in self.cogs_per_contig.df.columns, phylum.scgs)
        total = len(phylum.scgs)
        assert scgs
        rows = {}
        for cluster, group_df in self.cogs_per_contig.df.groupby('cluster'):
            summed_counts = group_df[scgs].sum()
            nonzero = summed_counts.count()
            exact_1 = summed_counts[summed_counts == 1].count()
            above_1 = nonzero - exact_1
            row = {"{0} total".format(phylum.name): total,
                   "{0} ==1".format(phylum.name): exact_1,
                   "{0} >1".format(phylum.name): above_1}
            rows[cluster] = row
        self.phyla.append(phylum)
        new_matrix = pd.DataFrame.from_dict(rows, orient='index')
        if self.result_matrix is None:
            self.result_matrix = new_matrix
        else:
            self.result_matrix = pd.concat([self.result_matrix, new_matrix], axis=1)
