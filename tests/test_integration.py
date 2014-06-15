#!/usr/bin/env python
from nose.tools import assert_equal, assert_true, assert_almost_equal, nottest
import nose.plugins.skip as skip
from os.path import isdir,isfile
from os import listdir
import os
import sys
import pandas as pd
from cogler import CogsPerContig, Phylum, Output
import subprocess

file_path = os.path.realpath(__file__)
test_dir_path = os.path.dirname(file_path)
test_data_path = os.path.join(test_dir_path, "test_data")
tmp_dir_path = os.path.join(test_dir_path, "nose_tmp_output")

class TestCogPhylumTable(object):
    """Test the COG_phylum_table script for cogler"""

    def setUp(self):
        """Create temporary dir if necessary,
        otherwise clear contents of it"""
        if not isdir(tmp_dir_path):
            os.mkdir(tmp_dir_path)
        self.tearDown()

    def tearDown(self):
        """remove temporary output files"""
        for f in os.listdir(tmp_dir_path):
            f_path = os.path.join(tmp_dir_path,f)
            try:
                os.remove(f_path)
            except:
                for e in os.listdir(f_path):
                    e_path = os.path.join(f_path,e)
                    os.remove(e_path)
                os.rmdir(f_path)
        assert os.listdir(tmp_dir_path) == []

    def test_correct_values_minimal(self):
        script_path = os.path.join(test_dir_path, "..", "scripts", "COG_phylum_table.py")
        blast_output = os.path.join(test_data_path, "minimal", "test_contigs.out")
        clustering_file = os.path.join(test_data_path, "minimal", "test_clustering_gt1000.csv")
        phylum_scg_file = os.path.join(test_data_path, "minimal", "test_phyla_scg.tsv")
        cdd_to_cog_file = os.path.join(test_dir_path, "..", "data", "cdd_to_cog.tsv")
        output_file = os.path.join(tmp_dir_path, "test_script_output_path.csv")
        
        cmd = ["python", script_path, 
                "-b", blast_output,
                "-c", clustering_file,
                "-m", phylum_scg_file,
                "--cdd_cog_file", cdd_to_cog_file,
                "--output_file", output_file]
        cmd = " ".join(cmd)
        output = subprocess.check_output(cmd, shell=True)
        
        result_matrix = pd.DataFrame.from_csv(output_file, index_col=0)
        assert_true(len(result_matrix[result_matrix > 0].columns) == 3)
        assert_true(result_matrix.ix[47, 'Acidobacteria total'] == 2)
        assert_true(result_matrix.ix[47, 'Acidobacteria >1'] == 1)
        assert_true(result_matrix.ix[47, 'Acidobacteria ==1'] == 1)

    def test_correct_values_medium(self):
        reason = "Test data not in place"
        raise skip.SkipTest(reason)

        script_path = os.path.join(test_dir_path, "..", "scripts", "COG_phylum_table.py")
        blast_output = os.path.join(test_data_path, "medium", "test_contigs.out")
        clustering_file = os.path.join(test_data_path, "medium", "test_clustering_gt1000.csv")
        phylum_scg_file = os.path.join(test_data_path, "medium", "test_phyla_scg.tsv")
        cdd_to_cog_file = os.path.join(test_dir_path, "..", "data", "cdd_to_cog.tsv")
        output_file = os.path.join(tmp_dir_path, "test_script_output_path.csv")
        
        cmd = ["python", script_path, 
                "-b", blast_output,
                "-c", clustering_file,
                "-m", phylum_scg_file,
                "--cdd_cog_file", cdd_to_cog_file,
                "--output_file", output_file]
        cmd = " ".join(cmd)
        output = subprocess.check_output(cmd, shell=True)
        
        result_matrix = pd.DataFrame.from_csv(output_file, index_col=0)

    def test_running_complex(self):
        reason = "Test data not in place"
        raise skip.SkipTest(reason)
        script_path = os.path.join(test_dir_path, "..", "scripts", "COG_phylum_table.py")
        blast_output = os.path.join(test_data_path, "test_contigs.out")
        clustering_file = os.path.join(test_data_path, "test_clustering_gt1000.csv")
        phylum_scg_file = os.path.join(test_data_path, "test_phyla_scg.tsv")
        cdd_to_cog_file = os.path.join(test_dir_path, "..", "data", "cdd_to_cog.tsv")
        output_file = os.path.join(tmp_dir_path, "test_script_output_path.csv")
        
        cmd = ["python", script_path, 
                "-b", blast_output,
                "-c", clustering_file,
                "-m", phylum_scg_file,
                "--cdd_cog_file", cdd_to_cog_file,
                "--output_file", output_file]
        cmd = " ".join(cmd)
        output = subprocess.check_output(cmd, shell=True)
        
        result_matrix = pd.DataFrame.from_csv(output_file, index_col=0)
        pass

