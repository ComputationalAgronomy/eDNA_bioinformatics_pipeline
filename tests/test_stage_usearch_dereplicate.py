import pytest
import os
import sys

from edna_processor.utils.base_logger import logger
from stage.stage_config import StageConfig
from stage.stage_usearch_dereplicate import DereplicateStage


@pytest.fixture
def config():
    config = StageConfig(verbose=True, dry=True, logger=logger)
    return config


@pytest.fixture
def stage(config):
    fasta_dir = "data_dir"
    save_dir = "output_dir"
    stage = DereplicateStage(config, fasta_dir=fasta_dir, save_dir=save_dir)
    return stage


def test_stage(stage):
    expected = "data_dir"
    assert stage.fasta_dir == expected

    expected = "output_dir"
    assert stage.save_dir == expected

    expected = "-sizeout -relabel Uniq"
    assert stage.params == expected


def test_setup(stage):
    stage.setup("test")

    summary = stage.summary()
    expected = ['Step 0: ==LOG== Program: Dereplicate trimmed sequences.', 'Step 1: ==LOG== RedirectOutput: Write usearch report.']
    assert summary == expected

    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test_cut.fasta")
    outfile = os.path.join("output_dir", "test_derep.fasta")
    expected = f"usearch.exe -fastx_uniques {infile} -sizeout -relabel Uniq -threads 1 -fastaout {outfile}"
    assert command == expected


def test_run(stage):
    dryrun = stage.run()
    output = """==LOG== Running: dereplicate\
    ==LOG== Program: usearch.exe."""
    assert dryrun == True


def test_params(config):
    fasta_dir = "data_dir"
    save_dir = "output_dir"
    annot_size = False
    seq_label = "Uniq_Test"

    stage = DereplicateStage(config, fasta_dir=fasta_dir, save_dir=save_dir,
                             annot_size=annot_size, seq_label=seq_label)

    expected = " -relabel Uniq_Test"
    assert stage.params == expected

    stage.setup("test2")
    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join(fasta_dir, "test2_cut.fasta")
    outfile = os.path.join(save_dir, "test2_derep.fasta")
    expected = f"usearch.exe -fastx_uniques {infile}  -relabel Uniq_Test -threads 1 -fastaout {outfile}"
    assert command == expected