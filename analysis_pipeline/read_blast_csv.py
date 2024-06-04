def read_hap_level(blast_csv_path):
    hap2level = {}

    print(f"> Blast CSV Table:  {blast_csv_path}")
    n = 0

    desired_level = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
    error_symbol= str.maketrans({':': '_', '/': '_', '\\': '_', '*': '_', '?': '_', '"': '_', '<': '_', '>': '_', '|': '_'})
    with open(blast_csv_path, 'r') as file:
        for line in file.readlines():
        
            line_list = line.split(',')
            haplotype = line_list[0]
            level_list = [str(line_list[i]) for i in range(2, 9)]
            # identity = line_list[9]
            # length = line_list[14]
            # evalue = line_list[17]
            # bitscore = line_list[18]

            level_list[0] = level_list[0].translate(error_symbol)
            if level_list[2] == 'Mugil':
                level_list[2] = "Mugilidae"

            hap2level[haplotype] = {desired_level[i]: level_list[i] for i in range(len(desired_level))}
            n += 1

    print(f"Haplotype Assigned:  {n} species\n")

    return hap2level