# eDNA_Bioinformatics_Pipeline
The purpose of this project is to implement eDNA bioinformatics processing and to explore taxonomic delimitation from various perspectives using eDNA sequences.

## Contents

- [Getting Started](#getting-started)
  - [Docker Version](#docker-version)
  - [Local Version](#local-version)
- [Pytest](#pytest)
- [Usage](#usage)
  - [`fastq_processor`](#fastq_processor-module)
  - [`analysis_toolkit`](#analysis_toolkit-module)

## Getting Started

### Docker Version

Download [Docker Desktop](https://docs.docker.com/get-docker/) on your own machine.

After installation, run `docker version` in terminal. If the version is displayed, it means the installation was successful.

Clone the GitHub repo:
```sh
git clone https://github.com/ComputationalAgronomy/eDNA_bioinformatics_pipeline.git
```

Build an image according to the Dockerfile.
```sh
docker build -t [ImageName] .
```

After the Dockerfile is successfully exported to an image, the [ImageName] should appear in the repo list if you use `docker image ls` to check.

Next, launch a new container from the Docker image that was just built:
```sh
docker run -it [ImageName]
```

If the launch is successful, your terminal should display something like the following, and **The Container is Ready to Work!**
```sh
(base) root@93f4d3cf355f:/#
```
`(base)` indicates the conda environment you are using, where all dependent Python packages are installed (don't deactivate this!). `93f4d3cf355f` indicates the ID of this container.

File structure inside a new container:
```
├── usr/
|   ├── local/bin/
|   |   ├── bbmap/
|   |   |   └── reformat.sh
|   |   ├── clustalo
|   |   ├── cutadapt
|   |   ├── iqtree-2.3.5-Linux-intel/
|   |   |   ├── bin/
|   |   |   |   └── iqtree2
|   |   ├── ncbi-blast-2.6.1+/
|   |   |   └── bin/
|   |   |       ├── blastn
|   |   |       └── makeblastdb
|   |   └── usearch
|   └── src/app/
|       ├── README.md
|       ├── dockerfile
|       ├── requirements.txt
|       ├── src/
|       |   ├── analysis_toolkit/
|       |   └── fastq_processor/
|       └── tests/
└── workplace/ # <- By default, you will start here.
    └── example/
```

#### Other useful commands when you are working with Docker:

(base) root@93f4d3cf355f:/# `exit`: Exit the container.

`docker cp [container ID]:/path/to/file /host/destination/folder`: Copying files from Docker container to host.

`docker cp /path/to/file [container ID]:/container/destination/folder`: Copying files from host to Docker container.

`docker container ls -a`: List all containers

`docker stop [container ID or Name]`: Stop a running container.

`docker start [container ID or Name]`: Start a stopped container.

`docker exec -it [container ID or Name] bash`: Enter a running container's shell.

`docker rmi [ImageName]`: Remove a Docker image.

`docker container rm [container ID or Name]`: Remove a container.

`docker system prune (--force)`: Remove \<none> TAG images (be careful when using this command).

### Local version

#### Dependency Installation

Make sure you have installed all of the following prerequisites on your machine:
* BBMap - [Download](https://sourceforge.net/projects/bbmap/) (You may also need to install java if you haven't already)
* Clustal Omega - [Download](http://www.clustal.org/omega/)
* Cutadapt - [Download](https://cutadapt.readthedocs.io/en/stable/installation.html)
* IQTREE2 - [Download](http://www.iqtree.org/)
* NCBI-BLAST+ - [Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* USEARCH - [Download](https://www.drive5.com/usearch/download.html)

, and ensure the path of downloaded software is added to the "Path" variable in "Environmental Variables".

#### Installation
1. Required Python Package Installation
```sh
python -m pip install -r requirements.txt
```
2. Local package installs
```sh
python -m pip install -e .
```

## Pytest
check `pytest.ini`
```sh
pytest
# OR
pytest tests
```

## Usage

### `fastq_processor` module
This module processes raw FASTQ data through several stages including paired-end merging, primer cutting, reformatting, dereplication, denoising, and taxonomic assignment.

**Simplist example**
```python
from fastq_processor import FastqProcessor

FastqProcessor(
    stages_parent_dir="./fastq_processor/stages",
    fastq_dir_name="fastq",
    db_path="./fastq_processor/db/MiFish",
    lineage_path="./fastq_processor/lineage/lineage.csv",
)
```

Ensure that `stages_parent_dir` exists and contains a subdirectory named `fastq_dir_name` with your raw FASTQ files.
```
stages/
└── fastq/
    ├── sample1_R1.fastq
    └── sample1_R2.fastq
```

The `db_path` serves as the reference data for assigning denoised sequences to species. It should be set to the folder path containing the indexed files, with the prefix string added. The MiFish index files provided in the *example* folder are built using the [complete+partial mtDNA sequence file](https://mitofish.aori.u-tokyo.ac.jp/species/detail/download/?filename=download%2F/complete_partial_mitogenomes.zip) downloaded from the MiFish Pipeline.

If you want to use a custom FASTA file as the reference database, you need to create the index using the `makeblastdb` command from ncbi-blast+. Here is the command:
```sh
makeblastdb -in ref.fasta -dbtype nucl -out db_prefix
```

The `lineage_path` serves as a reference for labeled species at the taxonomic level above the species. It should be set to the path of the lineage file. The [lineage.csv](https://github.com/billzt/MiFish/blob/main/mifish/data/lineage.csv) provided in the *example* folder is downloaded from the MiFish GitHub repository. If you want to use another lineage file, please ensure it includes taxonomic information from the domain to genus level and maintains the same format.

#### Other optional parameters:

`enabled_stages`:

`n_cpu`: Number of CPU cores to be used for processing. Default is `1`.

`verbose`: Set to True to enable detailed logging output. Default is `True`.

Paired-end merging related:

`maxdiff`

`pctid`

Primer cutting related:

`rm_p_5`

`rm_p_3`

`error_rate`

`min_read_len`

`max_read_len`

Denoising related:

`minsize`

`alpha`

Taxonomic assignment related:

`evalue`

`qcov_hsp_perc`

`perc_identity`

`outfmt`

`specifiers`

### `Analysis_toolkit` module