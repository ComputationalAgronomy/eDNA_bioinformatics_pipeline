import os
import re

folder = './data/keelung/test/6_blast'
filename_list = os.listdir(folder)

# for filename in filename_list:
#     if filename.endswith('_1.fastq'):
#         new_filename = filename.replace('_1.fastq', '_R1.fastq')
#         os.rename(folder+filename, folder+new_filename)
#     elif filename.endswith('_2.fastq'):
#         new_filename = filename.replace('_2.fastq', '_R2.fastq')
        # os.rename(folder+filename, folder+new_filename)

for filename in filename_list:
    new_filename = filename.replace('_zotu', '')
    os.rename(f'{folder}/{filename}', f'{folder}/{new_filename}')