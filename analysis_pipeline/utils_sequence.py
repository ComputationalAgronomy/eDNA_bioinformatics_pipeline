from Bio import SeqIO
import numpy as np
import os
import subprocess
import tempfile

def derep_fasta(seq_path: str, uniq_path: str, relabel: str, threads: int = 12, sizeout: bool = False) -> None:
    """
    Dereplicate a FASTA file by removing duplicate sequences.

    :param seq_path: Path to the input FASTA file.
    :param uniq_path: Path to the output FASTA file with unique sequences.
    :param relabel: Prefix to add to the sequence labels.
    :param threads: Number of threads to use. Default is 12.
    :param sizeout: If True, size annotations will be added to the output sequence labels. Default is False.
    """
    cmd = [
        'usearch', '-fastx_uniques', seq_path, '-threads', str(threads),
        '-relabel', f'{relabel}-', '-fastaout', uniq_path
    ]
    if sizeout:
        cmd.append('-sizeout')

    print("> Running USEARCH command:", ' '.join(cmd))
    try:
        subprocess.run(cmd, check=True)
        print(f"\n> Dereplicated FASTA file saved to: {uniq_path}\n")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during dereplication: {e}")

def write_fasta(units2fasta_dict: dict[str, str], seq_path: str, dereplicate: bool = False, sizeout: bool = False) -> None:
    """
    Write sequences to a FASTA file.

    :param units2fasta_dict: Dictionary with sequence names as keys and FASTA format sequences as values.
    :param seq_path: Path to the output FASTA file.
    :param dereplicate: If True, dereplicate the sequences before writing to the file. Default is False.
    :param sizeout: If True, size annotations will be added to the output sequence labels. Only works when dereplicating. Default is False.
    """
    if dereplicate:
        fasta_list = []
        with tempfile.TemporaryDirectory() as tmpdirname:
            for unit_name, unit_seq in units2fasta_dict.items():
                unit_fa_path = os.path.join(tmpdirname, f'{unit_name}.fa')
                unit_uniq_fa_path = os.path.join(tmpdirname, f'{unit_name}_uniq.fa')

                with open(unit_fa_path, 'w') as file:
                    file.write(unit_seq)

                derep_fasta(seq_path=unit_fa_path, uniq_path=unit_uniq_fa_path, relabel=unit_name, sizeout=sizeout)

                with open(unit_uniq_fa_path, 'r') as file:
                    fasta_list.append(file.read())
    else:
        fasta_list = list(units2fasta_dict.values())

    fasta_str = "".join(fasta_list)
    with open(seq_path, 'w') as file:
        file.write(fasta_str)

    num_seq = fasta_str.count(">")
    print(f"\n> Written {num_seq} sequences to: {seq_path}")

    return num_seq

def align_fasta(seq_path: str, aln_path: str) -> None:
    """
    Align sequences in a FASTA file using Clustal Omega and save the aligned sequences to an output file.

    :param seq_path: Path to the input FASTA file containing sequences to align.
    :param aln_path: Path to the output FASTA file to save the aligned sequences.
    """
    cmd = [
        'clustalo', '-i', seq_path, '-o', aln_path, '--force'
    ]
    print("Running Clustal Omega command:", ' '.join(cmd))

    try:
        subprocess.run(cmd, check=True)
        print(f"\n> Aligned FASTA file saved to: {aln_path}\n")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during alignment: {e}")

def get_uniq_seq_freq(seqs_path: str, uniq_seqs_path: str, seq_labels: list[str]) -> str:
    """
    Count the frequency of each label category for each unique sequence in the 'uniq_seqs_path' file based on the 'seqs_path' file.

    :param seqs_path: Path to the input FASTA file containing all sequences.
    :param uniq_seqs_path: Path to the input FASTA file containing only unique sequences.
    :param seq_labels: List of labels corresponding to each sequence in the 'seqs_path' file.
    :return: String in the NEXUS format required by PopART to build a haplotype network.
    """
    uniq_labels = np.unique(seq_labels)
    label_freq_each_uniq_seq = {}

    with open(uniq_seqs_path, 'r') as uniq_handle:
        uniq_records = list(SeqIO.parse(uniq_handle, 'fasta'))
        for record in uniq_records:
            label_freq_each_uniq_seq[record.name] = {uniq_label: 0 for uniq_label in uniq_labels}

    with open(seqs_path, 'r') as seq_handle:
        for i, seq_record in enumerate(SeqIO.parse(seq_handle, 'fasta')):
            for uniq_record in uniq_records:
                if seq_record.seq == uniq_record.seq:
                    label_freq_each_uniq_seq[uniq_record.name][seq_labels[i]] += 1

    freq_string = (
        f"Begin Traits;\n"
        f"Dimensions NTraits={len(uniq_labels)};\n"
        f"Format labels=yes missing=? separator=Comma;\n"
        f"TraitLabels {' '.join(uniq_labels)};\n"
        f"Matrix\n"
    )

    for seq_id, label_freq in label_freq_each_uniq_seq.items():
        freq_values = ",".join(map(str, label_freq.values()))
        freq_string += f"{seq_id} {freq_values}\n"

    freq_string += ";\nend;\n"

    return freq_string