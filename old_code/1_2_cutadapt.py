import os

def cut_adapt(merge_dir, save_dir, cpu=2, rm_p_5='GTCGGTAAAACTCGTGCCAGC', rm_p_3='CAAACTGGGATTAGATACCCCACTATG', min_read_len=204, max_read_len=254):
    files = os.listdir(merge_dir)
    number = []
    for filename in files:
        if filename.endswith('_merge.fastq'):
            number.append(filename.replace('_merge.fastq', ''))
    for num in number: 
        cmd = (f'cutadapt -g "{rm_p_5};max_error_rate=0.15...{rm_p_3};max_error_rate=0.15" \
                {merge_dir}{num}_merge.fastq --discard-untrimmed -j {cpu} \
                -m {min_read_len-len(rm_p_5)-len(rm_p_3)} -M {max_read_len-len(rm_p_5)-len(rm_p_3)} \
                >{save_dir}{num}_cut.fastq 2>{save_dir}{num}_report.txt')
        os.system(cmd)
