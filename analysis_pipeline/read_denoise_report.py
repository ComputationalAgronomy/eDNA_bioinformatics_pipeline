import re

from analysis_pipeline.read_blast_csv import Reader

class ReadDenoiseReport(Reader):
    RE_DENOISE_PATTERN = re.compile(r"(Uniq\d*)|size=(\d*)|(amp\d*)|top=(Uniq\d*)")
    RE_CHFILTER_PATTERN = re.compile(r"(Uniq\d*)|size=(\d*)|(zotu)|(chimera)")

    def __init__(self):
        super().__init__()
        self.amp_size = {}
        self.hap2amp = {}
        self.hap_size = {}

# TODO(SW): Make this a instance method later
def read_denoise_report(zotu_report_path):
    amp_size = {}
    hap2amp = {}
    hap_size = {}

    print(f"> Denoising Report:  {zotu_report_path}")

    zotu_n = 0
    ch_n = 0

    with open(zotu_report_path, 'r') as file:
        for line in file.readlines():
            # create relationship between haplotypes(top) and amplicons
            if 'denoise' in line:
                # TODO(SW): If you refactor this into functions, you can make this more readable and easier to do pytest on it.
                line_list = [
                    "".join(t)
                    for t in ReadDenoiseReport.RE_DENOISE_PATTERN.findall(line)
                ]
                # e.g of line_list: ['Uniq1', '88422', 'amp1'], ['Uniq6', '8126', 'Uniq2'] <- [Amplicon_id, Size, Top]. If Top is "ampXX", it is a true amplicon, otherwise, it is a noise.
                if 'amp' in line_list[2]:
                    line_list[2] = line_list[0]

                amplicon, size, top = line_list[0], line_list[1], line_list[2]

                amp_size[amplicon] = size

                if top not in hap2amp:
                    hap2amp[top] = []
                hap2amp[top].append(amplicon)

            # rename top
            elif 'chfilter' in line:
                line_list = [
                    "".join(t)
                    for t in ReadDenoiseReport.RE_CHFILTER_PATTERN.findall(line)
                ]
                # e.g of line_list: ['Uniq1', '88422', 'zotu'], ['Uniq102', '124', 'chimera'] <- [Amplicon_id, assigned_type]
                old_top, size, assigned_type = line_list[0], line_list[1], line_list[2]

                if assigned_type == 'zotu':
                    zotu_n += 1 # TODO(SW): So you will be starting from 1 instead of 0?
                    new_top = f'Zotu{zotu_n}'
                    hap_size[new_top] = size
                elif assigned_type == 'chimera':
                    ch_n += 1
                    new_top = f'Chimera{ch_n}'

                hap2amp[new_top] = hap2amp.pop(old_top)

    print(f"Biological Haplotypes:  {zotu_n} kept\nPredicted Chimeras:  {ch_n} removed")

    return amp_size, hap2amp, hap_size
