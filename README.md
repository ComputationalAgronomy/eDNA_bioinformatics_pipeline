# eDNA_Bioinformatics_Pipeline
The purpose of this project is to implement eDNA bioinformatics processing and to explore taxonomic delimitation from various perspectives using eDNA sequences.

## Getting Started

### Docker Version
```sh
docker build -t [ImageName] .
```

After the Dockerfile is successfully exported to an image, the [ImageName] should appear in the repo list if you use `docker image ls` to check.

```sh
docker image ls
```
### Local version

#### Prerequisites
Make sure you have installed all of the following prerequisites on your machine:
* USEARCH - [Download](https://www.drive5.com/usearch/download.html)
* Cutadapt - [Download](https://cutadapt.readthedocs.io/en/stable/installation.html)
* NCBI-BLAST+ - [Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* Clustal Omega - [Download](http://www.clustal.org/omega/)
* IQTREE2 - [Download](http://www.iqtree.org/)

, and ensure the path of downloaded software is added to the "Path" variable in "Environmental Variables".

#### Installation
1. Clone the repo
```sh
git clone https://github.com/ComputationalAgronomy/eDNA_bioinformatics_pipeline.git
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

### Python Interface
```python
from import
```

### Command Line Interface
```sh

```