class Reader:
    def __init__(self):
        pass


class BlastReader(Reader):
    DESIRED_LEVEL = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

    @staticmethod
    def generate_error_table(error_code=':/\\*?"<>|', replace_symbol='_'):
        key2 = list(error_code)
        error_translation = dict.fromkeys(key2, replace_symbol)
        error_symbol = str.maketrans(error_translation)
        return error_symbol

    def __init__(self):
        super().__init__()
        self.error_table = BlastReader.generate_error_table()

    TAX_REPLACMENT = {"Mugil": "Mugilidae",
                      "KEY2": "REPLACE2"} # Allow multiple replacements later

# TODO(SW): Make this a instance method later
def read_hap_level(blast_csv_path):
    br = BlastReader()
    hap2level = {}

    print(f"> Blast CSV Table:  {blast_csv_path}")
    n = 0

    with open(blast_csv_path, 'r') as file:
        for line in file.readlines():
            line_list = line.split(',')
            haplotype = line_list[0]
            level_list = [str(line_list[i]) for i in range(2, 9)]
            # identity = line_list[9]
            # length = line_list[14]
            # evalue = line_list[17]
            # bitscore = line_list[18]

            level_list[0] = level_list[0].translate(br.error_table)
            if level_list[2] == BlastReader.TAX_REPLACMENT:
                level_list[2] = BlastReader.TAX_REPLACMENT[level_list[2]]

            hap2level[haplotype] = dict(zip(Reader.DESIRED_LEVEL, level_list))
            n += 1 # Will this be equal to the number of lines in teh file?

    print(f"Haplotype Assigned:  {n} species\n")

    return hap2level