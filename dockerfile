FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y wget gzip tar sudo python3 python3-pip python3-venv && \
    apt-get clean

WORKDIR /usr/local/bin

RUN wget https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz && \
    gunzip usearch11.0.667_i86linux32.gz && \
    chmod +x usearch11.0.667_i86linux32

RUN sudo apt-get update && \
    sudo apt-get install -y cutadapt

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz

RUN wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64 && \
    chmod +x clustalo-1.2.4-Ubuntu-x86_64

RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.3.5/iqtree-2.3.5-Linux-intel.tar.gz && \
    tar -zxvf iqtree-2.3.5-Linux-intel.tar.gz

WORKDIR /usr/src/app

COPY . .

RUN python3 -m venv .venv

RUN /bin/bash -c "source .venv/bin/activate && python3 -m pip install --upgrade pip && python3 -m pip install -r requirements.txt"

RUN /bin/bash -c "source .venv/bin/activate && python3 -m pip install -e ."

ENV PATH="/usr/local/bin/ncbi-blast-2.16.0+/bin:/usr/local/bin/iqtree-2.3.5-Linux-intel/bin:${PATH}"

CMD ["bash"]