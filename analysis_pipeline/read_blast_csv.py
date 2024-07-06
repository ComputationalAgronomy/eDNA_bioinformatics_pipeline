class Reader:
    def __init__(self):
        pass


class BlastReader(Reader):
    DESIRED_LEVEL = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]

    TAX_REPLACMENT = {"Mugil": "Mugilidae"}
                      #"KEY2": "REPLACE2"} # Allow multiple replacements later

    @staticmethod
    def generate_error_table(error_code=':/\\*?"<>|', replace_symbol='_'):
        """
        Generate an error table to translate illegal characters in species names to a standard symbol.

        :param error_code: The error code used in the BLAST table.
        :param replace_symbol: The symbol to replace the error code with.
        :return: A dictionary containing the error code as the key and the replace symbol as the value.
        """
        key2 = list(error_code)
        error_translation = dict.fromkeys(key2, replace_symbol)
        error_symbol = str.maketrans(error_translation)
        return error_symbol

    def __init__(self):
        super().__init__()
        self.error_table = BlastReader.generate_error_table()

    def read_hap_level(self, blast_csv_path: str) -> dict[str, dict[str, str]]:
        """
        Read the BLAST CSV table and return a dictionary of the corresponding taxonomic terms at each level for every haplotype (ZOTU).
        Seven levels are used: species, genus, family, order, class, phylum, kingdom.

        :param blast_csv_path: Path to the BLAST CSV table.
        :return: A dictionary containing the taxonomic terms for each haplotype. e.g. {"Zotu1": {"species": "spcA", ...}, ...}
        """
        print(f"> Blast CSV Table:  {blast_csv_path}")

        hap2level = {}
        with open(blast_csv_path, 'r') as file:
            for line in file.readlines():
                line_list = line.split(',')
                haplotype = line_list[0]
                level_list = [str(line_list[i]) for i in range(2, 9)]
                # identity = line_list[9]
                # length = line_list[14]
                # evalue = line_list[17]
                # bitscore = line_list[18]
    
                level_list[0] = level_list[0].translate(self.error_table)
                if level_list[2] in self.TAX_REPLACMENT:
                    level_list[2] = self.TAX_REPLACMENT[level_list[2]]

                hap2level[haplotype] = dict(zip(self.DESIRED_LEVEL, level_list))

        hap_count = len(hap2level)
        print(f"Haplotype Assigned:  {hap_count} species\n")
    
        return hap2level