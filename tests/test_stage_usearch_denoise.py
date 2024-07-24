import pytest
import os
import sys

from edna_processor.utils.base_logger import logger
from stage.stage_config import StageConfig
from stage.stage_usearch_denoise import DenoiseStage


@pytest.fixture
def config():
    config = StageConfig(verbose=True, dry=True, logger=logger)
    return config


@pytest.fixture
def stage(config):
    derep_dir = "data_dir"
    save_dir = "output_dir"
    stage = DenoiseStage(config, derep_dir=derep_dir, save_dir=save_dir)
    return stage


def test_stage(stage):
    expected = "data_dir"
    assert stage.derep_dir == expected

    expected = "output_dir"
    assert stage.save_dir == expected

    expected = "-minsize 8 -unoise_alpha 2"
    assert stage.params == expected


def test_setup(stage):
    stage.setup("test")

    summary = stage.summary()
    expected = ['Step 0: ==LOG== Program: Denoise unique sequences.', 'Step 1: ==LOG== RedirectOutput: Write usearch report.']
    assert summary == expected

    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test_derep.fasta")
    outfile = os.path.join("output_dir", "test_denoise.fasta")
    denoise_report = os.path.join("output_dir", "test_denoise_report.txt")
    expected = f"usearch.exe -unoise3 {infile} -minsize 8 -unoise_alpha 2 -threads 1 -zotus {outfile} -tabbedout {denoise_report}"
    assert command == expected


def test_run(stage):
    dryrun = stage.run()
    output = """==LOG== Running: denoise\
    ==LOG== Program: usearch.exe."""
    assert dryrun == True


def test_params(config):
    derep_dir = "data_dir"
    save_dir = "output_dir"
    minsize = 4
    alpha = 4

    stage = DenoiseStage(config, derep_dir=derep_dir, save_dir=save_dir,
                         minsize=minsize, alpha=alpha)

    expected = "-minsize 4 -unoise_alpha 4"
    assert stage.params == expected

    stage.setup("test2")
    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join(derep_dir, "test2_derep.fasta")
    outfile = os.path.join(save_dir, "test2_denoise.fasta")
    denoise_report = os.path.join(save_dir, "test2_denoise_report.txt")
    expected = f"usearch.exe -unoise3 {infile} -minsize 4 -unoise_alpha 4 -threads 1 -zotus {outfile} -tabbedout {denoise_report}"
    assert command == expected