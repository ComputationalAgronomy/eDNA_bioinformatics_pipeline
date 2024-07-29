import os
import pytest
import sys

from edna_processor.utils.base_logger import logger
from stage.stage_config import StageConfig
from stage.stage_cutadapt_cut_primer import CutPrimerStage



@pytest.fixture
def config():
    config = StageConfig(verbose=True, dry=True, logger=logger)
    return config


@pytest.fixture
def stage(config):
    merge_dir = "data_dir"
    save_dir = "output_dir"
    stage = CutPrimerStage(config, merge_dir=merge_dir, save_dir=save_dir)
    return stage


def test_stage(stage):
    expected = "data_dir"
    assert stage.merge_dir == expected

    expected = "output_dir"
    assert stage.save_dir == expected

    expected = "-g GTCGGTAAAACTCGTGCCAGC;max_error_rate=0.15...CAAACTGGGATTAGATACCCCACTATG;max_error_rate=0.15 --minimum-length 156 --maximum-length 206"
    assert stage.params == expected


def test_setup(stage):
    stage.setup("test")

    summary = stage.summary()
    expected = ['Step 0: ==LOG== Program: Cut primers for merged sequences.', 'Step 1: ==LOG== RedirectOutput: Write cutadapt output, report.']
                # 'Step 1: ==LOG== RedirectOutput: Redirect cutadapt output.']
    assert summary == expected

    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test_merge.fastq")
    expected = f"cutadapt {infile} -g GTCGGTAAAACTCGTGCCAGC;max_error_rate=0.15...CAAACTGGGATTAGATACCCCACTATG;max_error_rate=0.15 --minimum-length 156 --maximum-length 206 --discard-untrimmed -j 1"
    assert command == expected


def test_run(stage):
    dryrun = stage.run()
    output = """==LOG== Running: cut_primer
    ==LOG== Program: cutadapt."""
    assert dryrun == True


def test_params(config):
    merge_dir = "data_dir"
    save_dir = "output_dir"

    rm_p_5="AAA"
    rm_p_3="CCCC"
    error_rate=0.2
    min_read_len=100
    max_read_len=150

    stage = CutPrimerStage(config, merge_dir=merge_dir, save_dir=save_dir,
                          rm_p_5 = rm_p_5,
                          rm_p_3 = rm_p_3,
                          error_rate = error_rate,
                          min_read_len = min_read_len,
                          max_read_len = max_read_len)

    expected = "-g AAA;max_error_rate=0.2...CCCC;max_error_rate=0.2 --minimum-length 93 --maximum-length 143"
    assert stage.params == expected

    stage.setup("test2")
    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test2_merge.fastq")
    expected = f"cutadapt {infile} -g AAA;max_error_rate=0.2...CCCC;max_error_rate=0.2 --minimum-length 93 --maximum-length 143 --discard-untrimmed -j 1"
    assert command == expected
