#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
from os.path import isdir,isfile
from os import listdir
import os
import sys

from cogler.io import read_blast_output, read_markers_file, read_clustering_file

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
test_data_path = os.path.join(test_dir_path, "test_data")

class TestIO(object):
    """Test input output for cogler"""

    def test_read_per_phylum_scgs(self):
        summary, phyla_scgs = read_per_phylum_scgs(
                os.path.join(test_data_path, "test_phyla_scg.tsv"))
        # First row
        assert_true(summary.ix['Acidobacteria'].Number_genomes == 5)
        assert_true(summary.ix['Acidobacteria'].Number_genera == 5)
        assert_true(summary.ix['Acidobacteria'].Number_family == 2)
        assert_true(summary.ix['Acidobacteria'].Number_classes == 2)
        assert_true(summary.ix['Acidobacteria'].Number_SCG == 428)
       
        # Last row
        assert_true(summary.ix['Gammaproteobacteria'].Number_genomes == 276)
        assert_true(summary.ix['Gammaproteobacteria'].Number_genera == 76)
        assert_true(summary.ix['Gammaproteobacteria'].Number_family == 28)
        assert_true(summary.ix['Gammaproteobacteria'].Number_classes == 1)
        assert_true(summary.ix['Gammaproteobacteria'].Number_SCG == 143)
        
        assert_true(phyla_scgs.ix['Acidobacteria'].COG0001 == 0)
        assert_true(phyla_scgs.ix['Acidobacteria'].COG5665 == 0)
        assert_true(phyla_scgs.ix['Acidobacteria'].COG0002 == 1)
        
        assert_true(phyla_scgs.ix['Gammaproteobacteria'].COG0001 == 0)
        assert_true(phyla_scgs.ix['Gammaproteobacteria'].COG0012 == 1)
        assert_true(phyla_scgs.ix['Gammaproteobacteria'].COG5665 == 0)

        assert_true(phyla_scgs.ix['Acidobacteria'].sum() == summary.ix['Acidobacteria'].Number_SCG)
        assert_true(phyla_scgs.ix['Gammaproteobacteria'].sum() == summary.ix['Gammaproteobacteria'].Number_SCG)
    def test_blast(self):
        records, sseq_ids = read_blast_output(
                os.path.join(test_data_path, "test_contigs.out"))
        assert_true(len(records) > 0)
        assert_true(len(sseq_ids) > 0)

        # The test data contains 462 lines
        assert_true(len(records) == 462)
        # The first item in the blast output
        assert_true('224121' in sseq_ids)
        first_record = records[0]
        assert_true(first_record['qseqid'] == 'contig00001_1')
        assert_true(first_record['sseqid'] == 'gnl|CDD|224121')
        assert_true(first_record['pident'] == 39.34)
        assert_true(first_record['sstart'] == 10.0)
        assert_true(first_record['send'] == 68.0)
        assert_true(first_record['slen'] == 677.0)

        # The last item in the blast output
        first_record = records[-1]
        assert_true(first_record['qseqid'] == 'contig00009_55')
        assert_true(first_record['sseqid'] == 'gnl|CDD|223735')
        assert_true(first_record['pident'] == 46.78)
        assert_true(first_record['sstart'] == 3)
        assert_true(first_record['send'] == 173)
        assert_true(first_record['slen'] == 176)

    def test_markers_file(self):
        markers = read_markers_file(
                os.path.join(test_data_path, "test_markers_file.txt"))
        assert_true(len(markers) == 36)
        assert_true(markers[0] == 'COG0016')
        assert_true(markers[-1] == 'COG0130')

    def test_clustering_file(self):
        clusters, contigs_per_cluster = read_clustering_file(
                os.path.join(test_data_path, "test_clustering_gt1000.csv"))
        assert_true(len(clusters) == 2)
        assert_true('47' in clusters)
        assert_true('94' in clusters)
