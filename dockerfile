########## Ubuntu-based

FROM ubuntu:latest

########## Install utilized tools

RUN apt-get update && \
    apt-get install -y wget gzip tar sudo python3 python3-pip && \
    apt-get clean

########## Install prerequisites

WORKDIR /usr/local/bin

# USEARCH
RUN wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz && \
    gunzip usearch11.0.667_i86linux32.gz && \
    mv usearch11.0.667_i86linux32 usearch.exe && \
    chmod +x usearch.exe

# Cutadapt
RUN sudo apt-get update && \
    sudo apt-get install -y cutadapt

# NCBI-BLAST+
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz

# Clustal Omega
RUN wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64 && \
    mv clustalo-1.2.4-Ubuntu-x86_64 clustalo.exe && \
    chmod +x clustalo.exe

# IQTREE2
RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.3.5/iqtree-2.3.5-Linux-intel.tar.gz && \
    tar -zxvf iqtree-2.3.5-Linux-intel.tar.gz

ENV PATH="/usr/local/bin/ncbi-blast-2.16.0+/bin:/usr/local/bin/iqtree-2.3.5-Linux-intel/bin:${PATH}"

########## Install Miniconda

WORKDIR /tmp

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/opt/conda/bin:${PATH}"

RUN conda init

########## Install required python packages & edna_bp (eDNA bioinformatics pipeline) 

WORKDIR /usr/src/app

COPY . .

RUN /bin/bash -c "source activate && python3 -m pip install --upgrade pip && python3 -m pip install -r requirements.txt && python3 -m pip install -e ."

########## Build workplace

WORKDIR /

RUN mkdir workplace

WORKDIR /workplace

CMD ["bash"]