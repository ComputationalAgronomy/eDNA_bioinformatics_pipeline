import sys
import pytest

from stage.stage_config import StageConfig
from stage.usearch_merge_stage import UsearchMergeStage


@pytest.fixture
def config():
    config = StageConfig(verbose=True, dry=True, logger=sys.stdout)
    return config


@pytest.fixture
def stage(config):
    fastq_dir = "data_dir"
    save_dir = "output_dir"
    stage = UsearchMergeStage(config, fastq_dir=fastq_dir, save_dir=save_dir)
    return stage


def test_stage(stage):
    expected = "data_dir"
    assert stage.fastq_dir == expected

    expected = "output_dir"
    assert stage.save_dir == expected

    expected = " -fastq_maxdiffs 5 -fastq_pctid 90 -threads 1"
    assert stage.params == expected


def test_setup(stage):
    stage.setup("test")

    summary = stage.summary()
    expected = ["Step 0: ==LOG== Program: usearch_merge."]
    assert summary == expected

    runner = stage.runners[0]
    command = runner.command
    expected = "usearch  -fastq_mergepairs data_dir/test_R1.fastq -fastqout data_dir/test_merged.fastq -fastq_maxdiffs 5 -fastq_pctid 90 -threads 1 -report output_dir/test_report.txt"
    assert command == expected


def test_run(stage):
    dryrun = stage.run()
    output = """==LOG== Running: usearch
    ==LOG== Program: usearch_merge."""
    assert dryrun == True


def test_params(config):
    fastq_dir = "data_dir"
    save_dir = "output_dir"
    maxdiff = 10
    pctid = 80

    stage = UsearchMergeStage(config, fastq_dir=fastq_dir, save_dir=save_dir,
                              maxdiff=maxdiff, pctid=pctid)

    expected = " -fastq_maxdiffs 10 -fastq_pctid 80 -threads 1"
    assert stage.params == expected

    stage.setup("test2")
    runner = stage.runners[0]
    command = runner.command
    expected = "usearch  -fastq_mergepairs data_dir/test2_R1.fastq -fastqout data_dir/test2_merged.fastq -fastq_maxdiffs 10 -fastq_pctid 80 -threads 1 -report output_dir/test2_report.txt"
    assert command == expected
