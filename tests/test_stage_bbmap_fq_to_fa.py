import pytest
import os
import sys

from stage.stage_config import StageConfig
from stage.stage_bbmap_fq_to_fa import FqToFaStage


@pytest.fixture
def config():
    config = StageConfig(verbose=True, dry=True, logger=sys.stdout)
    return config


@pytest.fixture
def stage(config):
    cutprimer_dir = "data_dir"
    save_dir = "output_dir"
    stage = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=save_dir)
    return stage


def test_stage(stage):
    expected = "data_dir"
    assert stage.cutprimer_dir == expected

    expected = "output_dir"
    assert stage.save_dir == expected

    expected = "overwrite=True"
    assert stage.params == expected


def test_setup(stage):
    stage.setup("test")

    summary = stage.summary()
    expected = ["Step 0: ==LOG== Program: bbmap_reformat.sh."]
    assert summary == expected

    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test_cut.fastq")
    outfile = os.path.join("output_dir", "test_cut.fasta")
    report = os.path.join("output_dir", "test_report.txt")
    expected = f"bash reformat.sh in={infile} out={outfile} overwrite=True 2>{report}"
    assert command == expected


def test_run(stage):
    dryrun = stage.run()
    output = """==LOG== Running: fq_to_fa\
    ==LOG== Program: bbmap_reformat.sh."""
    assert dryrun == True


def test_params(config):
    cutprimer_dir = "data_dir"
    save_dir = "output_dir"
    overwrite = False

    stage = FqToFaStage(config, cutprimer_dir=cutprimer_dir, save_dir=save_dir,
                        overwrite=overwrite)

    expected = "overwrite=False"
    assert stage.params == expected

    stage.setup("test2")
    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test2_cut.fastq")
    outfile = os.path.join("output_dir", "test2_cut.fasta")
    report = os.path.join("output_dir", "test2_report.txt")
    expected = f"bash reformat.sh in={infile} out={outfile} overwrite=False 2>{report}"
    assert command == expected
