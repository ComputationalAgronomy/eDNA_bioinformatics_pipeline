import pytest
from analysis_pipeline.read_blast_csv import BlastReader

class TestBlastReader:

    def test_generate_error_table(self):
        reader = BlastReader()
        error_symbol = reader.generate_error_table()
        expected_output = {58: '_', 47: '_', 92: '_', 42: '_', 63: '_', 34: '_', 60: '_', 62: '_', 124: '_'}

        assert error_symbol == expected_output
        assert reader.error_table == expected_output

    def test_generate_error_table(self):
        reader = BlastReader()
        error_symbol = reader.generate_error_table(error_code='abcABC', replace_symbol='=')
        expected_output = {97: '=', 98: '=', 99: '=', 65: '=', 66: '=', 67: '='}

        assert error_symbol == expected_output
