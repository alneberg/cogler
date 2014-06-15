#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
import nose.plugins.skip as skip
from os.path import isdir,isfile
from os import listdir
import os
import sys
import pandas as pd
from cogler import CogsPerContig, Phylum, Output

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
test_data_path = os.path.join(test_dir_path, "test_data")

class TestCogsPerClusters(object):
    """Test the CogsPerClusters class for cogler"""
    
    def setUp(self):
        self.cluster_per_contig = {'contig1': '1',
                                   'contig2': '2',
                                   'contig3': '1'}

        self.features_per_contig = {'contig1': ['COG0001', 'COG0004'],                           
                                    'contig2': ['COG0002', 'COG0004'],
                                    'contig3': ['COG0003', 'COG0001']} 

        self.cogs_per_contig = CogsPerContig(self.cluster_per_contig,
                                             self.features_per_contig)
    
    def test_init(self):
        """Test for the CogsPerContig init method."""
        assert_true(self.cogs_per_contig.cluster_per_contig == self.cluster_per_contig)
        assert_true(self.cogs_per_contig.features_per_contig == self.features_per_contig)

    def test_correct_values(self):
        """Test that the values in CogPerContig are correct."""
        assert_true(self.cogs_per_contig.df.ix['contig1','COG0001'] == 1)
        assert_true(self.cogs_per_contig.df.ix['contig1','COG0004'] == 1)
        assert_true(self.cogs_per_contig.df.ix['contig2','COG0002'] == 1)
        assert_true(self.cogs_per_contig.df.ix['contig2','COG0004'] == 1)
        assert_true(self.cogs_per_contig.df.ix['contig3','COG0003'] == 1)
        assert_true(self.cogs_per_contig.df.ix['contig3','COG0001'] == 1)

        # The nulls
        assert_true(pd.isnull(self.cogs_per_contig.df.ix['contig1','COG0002']))
        assert_true(pd.isnull(self.cogs_per_contig.df.ix['contig1','COG0003']))
        assert_true(pd.isnull(self.cogs_per_contig.df.ix['contig2','COG0001']))
        assert_true(pd.isnull(self.cogs_per_contig.df.ix['contig2','COG0003']))
        assert_true(pd.isnull(self.cogs_per_contig.df.ix['contig3','COG0002']))
        assert_true(pd.isnull(self.cogs_per_contig.df.ix['contig3','COG0004']))

    def test_cluster_membership(self):
        """Assert the contigs have the correct cluster assignment."""
        assert_true(self.cogs_per_contig.df.ix['contig1','cluster'] == '1')
        assert_true(self.cogs_per_contig.df.ix['contig2','cluster'] == '2')
        assert_true(self.cogs_per_contig.df.ix['contig3','cluster'] == '1')


    def test_number_of_values(self):
        """Assert nothing else is in the dataframe."""
        assert_true(self.cogs_per_contig.df.shape == (3,5))


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
        self.cluster_per_contig = {'contig1': '1',
                                   'contig2': '2',
                                   'contig3': '1'}

        self.features_per_contig = {'contig1': ['COG0001', 'COG0004'],
                                    'contig2': ['COG0002', 'COG0004'],
                                    'contig3': ['COG0003', 'COG0001']}
        self.cogs_per_contig = CogsPerContig(self.cluster_per_contig,
                                                 self.features_per_contig)

        self.phylum1 = Phylum("Phylum1", ['COG0001', 'COG0004'])
        self.phylum2 = Phylum("Phylum2", ['COG0002', 'COG0003', 'COG0004'])
        self.phylum3 = Phylum("Phylum3", ['COG0010', 'COG0004'])

    def test_simple(self):
        """A simple use case test for the Output class"""
        output = Output(self.cogs_per_contig)
        output.add_phylum(self.phylum1)
        assert_true(set(output.result_matrix.columns) == set(['Phylum1 total', 'Phylum1 >1', 'Phylum1 ==1']))
        assert_true(len(output.result_matrix.index.values) == 2)
        assert_true(set(output.result_matrix.index.values) == set(['1', '2']))
        assert_true(output.result_matrix.ix['1', 'Phylum1 total'] == 2)
        assert_true(output.result_matrix.ix['1', 'Phylum1 >1'] == 1)
        assert_true(output.result_matrix.ix['1', 'Phylum1 ==1'] == 1)

        assert_true(output.result_matrix.ix['2', 'Phylum1 total'] == 2)
        assert_true(output.result_matrix.ix['2', 'Phylum1 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum1 ==1'] == 1)

    def test_all(self):
        """A more complex use case with three phyla for the Output class"""
        output = Output(self.cogs_per_contig)
        output.add_phylum(self.phylum1)
        output.add_phylum(self.phylum2)
        output.add_phylum(self.phylum3)

        columns = ['Phylum1 total', 'Phylum1 >1', 'Phylum1 ==1',
                'Phylum2 total', 'Phylum2 >1', 'Phylum2 ==1',
                'Phylum3 total', 'Phylum3 >1', 'Phylum3 ==1']

        for column in columns:
            assert_true(column in output.result_matrix.columns)
        assert_true(len(output.result_matrix.index.values) == 2)
        assert_true(set(output.result_matrix.index.values) == set(['1', '2']))
        assert_true(output.result_matrix.ix['1', 'Phylum1 total'] == 2)
        assert_true(output.result_matrix.ix['1', 'Phylum1 >1'] == 1)
        assert_true(output.result_matrix.ix['1', 'Phylum1 ==1'] == 1)

        assert_true(output.result_matrix.ix['2', 'Phylum1 total'] == 2)
        assert_true(output.result_matrix.ix['2', 'Phylum1 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum1 ==1'] == 1)

        assert_true(output.result_matrix.ix['1', 'Phylum2 total'] == 3)
        assert_true(output.result_matrix.ix['1', 'Phylum2 >1'] == 0)
        assert_true(output.result_matrix.ix['1', 'Phylum2 ==1'] == 2)

        assert_true(output.result_matrix.ix['2', 'Phylum2 total'] == 3)
        assert_true(output.result_matrix.ix['2', 'Phylum2 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum2 ==1'] == 2)

        assert_true(output.result_matrix.ix['1', 'Phylum3 total'] == 2)
        assert_true(output.result_matrix.ix['1', 'Phylum3 >1'] == 0)
        assert_true(output.result_matrix.ix['1', 'Phylum3 ==1'] == 1)
        
        assert_true(output.result_matrix.ix['2', 'Phylum3 total'] == 2)
        assert_true(output.result_matrix.ix['2', 'Phylum3 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum3 ==1'] == 1)

    def test_multiple_on_contig(self):
        """The output should not be affected if the same cog is present 
        twice in one contig compared to only once"""
        features_per_contig = self.features_per_contig
        features_per_contig['contig3'].append('COG0003')
        self.cogs_per_contig = CogsPerContig(self.cluster_per_contig,
                features_per_contig)
        output = Output(self.cogs_per_contig)
        output.add_phylum(self.phylum1)
        output.add_phylum(self.phylum2)
        output.add_phylum(self.phylum3)

        columns = ['Phylum1 total', 'Phylum1 >1', 'Phylum1 ==1',
                'Phylum2 total', 'Phylum2 >1', 'Phylum2 ==1',
                'Phylum3 total', 'Phylum3 >1', 'Phylum3 ==1']

        for column in columns:
            assert_true(column in output.result_matrix.columns)
        assert_true(len(output.result_matrix.index.values) == 2)
        assert_true(set(output.result_matrix.index.values) == set(['1', '2']))
        assert_true(output.result_matrix.ix['1', 'Phylum1 total'] == 2)
        assert_true(output.result_matrix.ix['1', 'Phylum1 >1'] == 1)
        assert_true(output.result_matrix.ix['1', 'Phylum1 ==1'] == 1)

        assert_true(output.result_matrix.ix['2', 'Phylum1 total'] == 2)
        assert_true(output.result_matrix.ix['2', 'Phylum1 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum1 ==1'] == 1)

        assert_true(output.result_matrix.ix['1', 'Phylum2 total'] == 3)
        assert_true(output.result_matrix.ix['1', 'Phylum2 >1'] == 0)
        assert_true(output.result_matrix.ix['1', 'Phylum2 ==1'] == 2)

        assert_true(output.result_matrix.ix['2', 'Phylum2 total'] == 3)
        assert_true(output.result_matrix.ix['2', 'Phylum2 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum2 ==1'] == 2)

        assert_true(output.result_matrix.ix['1', 'Phylum3 total'] == 2)
        assert_true(output.result_matrix.ix['1', 'Phylum3 >1'] == 0)
        assert_true(output.result_matrix.ix['1', 'Phylum3 ==1'] == 1)
        
        assert_true(output.result_matrix.ix['2', 'Phylum3 total'] == 2)
        assert_true(output.result_matrix.ix['2', 'Phylum3 >1'] == 0)
        assert_true(output.result_matrix.ix['2', 'Phylum3 ==1'] == 1)

