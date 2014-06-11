#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
import nose.plugins.skip as skip
from os.path import isdir,isfile
from os import listdir
import os
import sys

from cogler import CogsPerClusters, Phylum, Output

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
test_data_path = os.path.join(test_dir_path, "test_data")

class TestCogsPerClusters(object):
    """Test the CogsPerClusters class for cogler"""
    
    def setUp(self):
        self.contigs_per_cluster = {'1': ['contig1', 'contig3'],
                '2': ['contig2']}
        self.features_per_contig = {'contig1': ['COG0001', 'COG0004'],                           
                                    'contig2': ['COG0002', 'COG0004'],
                                    'contig3': ['COG0003', 'COG0001']} 

        self.cogs_per_clusters = CogsPerClusters(self.contigs_per_clusters, 
                                                self.features_per_contig)
    
    def test_init(self):
        """Test for the CogPerCluster init method."""
        assert_true(self.cogs_per_clusters.contigs_per_cluster == self.contigs_per_cluster)
        assert_true(self.cogs_per_clusters.features_per_contig == self.features_per_contig)

    def test_count_features_per_cluster(self):
        """Test for count_features_per_cluster method."""
        features_per_cluster = self.cogs_per_clusters.count_features_per_cluster()

        assert_true(features_per_cluster['1'] == {'COG0001': 2, 'COG0003': 1, 'COG0004': 1})
        assert_true(features_per_cluster['2'] == {'COG0002': 1, 'COG0004': '1'})

        assert_true(len(features_per_clusters) == 2)

class TestPhylum(object):
    """Test the Phylum class for cogler"""
    def setUp(self):
        self.name = "Phylum"
        self.scgs = ['COG0001', 'COG0004']
        self.phylum = Phylum(self.name, self.scgs)

    def test_init(self):
        assert_true(self.phylum.name == self.name)
        assert_true(self.phylum.scgs == self.scgs)

class TestOutput(object):
    """Test the Output class for cogler"""
    def setUp(self):
        self.contigs_per_cluster = {'1': ['contig1', 'contig3'],                        
                '2': ['contig2']}

        self.features_per_contig = {'contig1': ['COG0001', 'COG0004'],
                                    'contig2': ['COG0002', 'COG0004'],
                                    'contig3': ['COG0003', 'COG0001']}
        self.cogs_per_clusters = CogsPerClusters(self.contigs_per_clusters,
                                                 self.features_per_contig)

        self.phylum1 = Phylum("Phylum1", ['COG0001', 'COG00004'])
        self.phylum2 = Phylum("Phylum2", ['COG0002', 'COG00003', 'COG0004'])
        self.phylum3 = Phylum("Phylum3", ['COG0010', 'COG0004'])

    def test_simple(self):
        """A simple use case test for the Output class"""
        output = Output(self.cogs_per_clusters)
        output.add_phylum(self.phylum1)
        assert_true(output.result_matrix.columns == ['cluster', 'Phylum1 total', 'Phylum1 >1', 'Phylum1 ==1'])
        assert_true(output.result_matrix.cluster == ['1', '2'])
        assert_true(output.result_matrix.ix[['1', 'Phylum1 total']] == 2)
        assert_true(output.result_matrix.ix[['1', 'Phylum1 >1']] == 1)
        assert_true(output.result_matrix.ix[['1', 'Phylum1 ==1']] == 1)

        assert_true(output.result_matrix.ix[['2', 'Phylum1 total']] == 2)
        assert_true(output.result_matrix.ix[['2', 'Phylum1 >1']] == 0)
        assert_true(output.result_matrix.ix[['2', 'Phylum1 ==1']] == 1)

    def test_all(self):
        """A more complex use case with three phyla for the Output class"""
        output = Output(self.cogs_per_clusters)
        output.add_phylum(self.phylum1)
        output.add_phylum(self.phylum2)
        output.add_phylum(self.phylum3)

        columns = ['cluster', 'Phylum1 total', 'Phylum1 >1', 'Phylum1 ==1',
                'Phylum2 total', 'Phylum2 >1', 'Phylum2 ==1',
                'Phylum3 total', 'Phylum3 >1', 'Phylum3 ==1']

        for column in columns:
            assert_true(column in output.result_matrix.columns)

        assert_true(output.result_matrix.cluster == ['1', '2'])
        assert_true(output.result_matrix.ix[['1', 'Phylum1 total']] == 2)
        assert_true(output.result_matrix.ix[['1', 'Phylum1 >1']] == 1)
        assert_true(output.result_matrix.ix[['1', 'Phylum1 ==1']] == 1)

        assert_true(output.result_matrix.ix[['2', 'Phylum1 total']] == 2)
        assert_true(output.result_matrix.ix[['2', 'Phylum1 >1']] == 0)
        assert_true(output.result_matrix.ix[['2', 'Phylum1 ==1']] == 1)

        assert_true(output.result_matrix.ix[['1', 'Phylum2 total']] == 3)
        assert_true(output.result_matrix.ix[['1', 'Phylum2 >1']] == 0)
        assert_true(output.result_matrix.ix[['1', 'Phylum2 ==1']] == 2)

        assert_true(output.result_matrix.ix[['2', 'Phylum2 total']] == 3)
        assert_true(output.result_matrix.ix[['2', 'Phylum2 >1']] == 0)
        assert_true(output.result_matrix.ix[['2', 'Phylum2 ==1']] == 2)

        assert_true(output.result_matrix.ix[['1', 'Phylum3 total']] == 2)
        assert_true(output.result_matrix.ix[['1', 'Phylum3 >1']] == 0)
        assert_true(output.result_matrix.ix[['1', 'Phylum3 ==1']] == 1)
        
        assert_true(output.result_matrix.ix[['2', 'Phylum3 total']] == 2)
        assert_true(output.result_matrix.ix[['2', 'Phylum3 >1']] == 0)
        assert_true(output.result_matrix.ix[['2', 'Phylum3 ==1']] == 1)


