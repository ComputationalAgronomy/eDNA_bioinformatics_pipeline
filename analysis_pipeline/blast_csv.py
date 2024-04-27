def read_hap_rank(blast_csv_path):
    hap2rank = {}

    print(f"> Blast CSV Table:  {blast_csv_path}")
    n = 0

    desired_rank = ["species", "genus", "family", "order", "class", "phylum", "kingdom"]
    error_symbol= str.maketrans({':': '_', '/': '_', '\\': '_', '*': '_', '?': '_', '"': '_', '<': '_', '>': '_', '|': '_'})
    with open(blast_csv_path, 'r') as file:
        for line in file.readlines():
        
            line_list = line.split(',')
            haplotype = line_list[0]
            rank_list = [str(line_list[i]) for i in range(2, 9)]
            # identity = line_list[9]
            # length = line_list[14]
            # evalue = line_list[17]
            # bitscore = line_list[18]

            rank_list[0] = rank_list[0].translate(error_symbol)
            if rank_list[2] == 'Mugil':
                rank_list[2] = "Mugilidae"

            hap2rank[haplotype] = {desired_rank[i]: rank_list[i] for i in range(len(desired_rank))}
            n += 1

    print(f"Haplotype Assigned:  {n} species\n")

    return hap2rank