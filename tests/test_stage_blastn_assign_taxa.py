import pytest
import os
import sys

from stage.stage_config import StageConfig
from stage.stage_blastn_assign_taxa import AssignTaxaStage


@pytest.fixture
def config():
    config = StageConfig(verbose=True, dry=True, logger=sys.stdout)
    return config


@pytest.fixture
def stage(config):
    denoise_dir = "data_dir"
    save_dir = "output_dir"
    db_path = "db_path"
    stage = AssignTaxaStage(config, denoise_dir=denoise_dir, save_dir=save_dir, db_path=db_path)
    return stage


def test_stage(stage):
    expected = "data_dir"
    assert stage.denoise_dir == expected

    expected = "output_dir"
    assert stage.save_dir == expected

    expected = f"-db db_path -max_target_seqs 1 -evalue 1e-05 -qcov_hsp_perc 90 -perc_identity 90 -outfmt \"10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\""
    assert stage.params == expected


def test_setup(stage):
    stage.setup("test")

    summary = stage.summary()
    expected = ["Step 0: ==LOG== Program: ncbi-blast+_blastn."]
    assert summary == expected

    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join("data_dir", "test_denoise.fasta")
    outfile = os.path.join("output_dir", "test_blast.csv")
    expected = f"blastn -query {infile} -db db_path -max_target_seqs 1 -evalue 1e-05 -qcov_hsp_perc 90 -perc_identity 90 -outfmt \"10 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -num_threads 1 -out {outfile}"
    assert command == expected


def test_run(stage):
    dryrun = stage.run()
    output = """==LOG== Running: denoise\
    ==LOG== Program: usearch.exe."""
    assert dryrun == True


def test_params(config):
    denoise_dir = "data_dir"
    save_dir = "output_dir"
    db_path = "db_path_test"
    maxhitnum = 5
    evalue = 1e-30
    qcov_hsp_perc = 85
    perc_identity = 97
    outfmt = "6"
    specifiers = "qseqid staxids"

    stage = AssignTaxaStage(config, denoise_dir=denoise_dir, save_dir=save_dir,
                            db_path=db_path, maxhitnum=maxhitnum, evalue=evalue,
                            qcov_hsp_perc=qcov_hsp_perc, perc_identity=perc_identity,
                            outfmt=outfmt, specifiers=specifiers)

    expected = "-db db_path_test -max_target_seqs 5 -evalue 1e-30 -qcov_hsp_perc 85 -perc_identity 97 -outfmt \"6 qseqid staxids\""
    assert stage.params == expected

    stage.setup("test2")
    runner = stage.runners[0]
    command = runner.command
    infile = os.path.join(denoise_dir, "test2_denoise.fasta")
    outfile = os.path.join(save_dir, "test2_blast.csv")
    expected = f"blastn -query {infile} -db db_path_test -max_target_seqs 5 -evalue 1e-30 -qcov_hsp_perc 85 -perc_identity 97 -outfmt \"6 qseqid staxids\" -num_threads 1 -out {outfile}"
    assert command == expected