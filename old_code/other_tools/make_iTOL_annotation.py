def read_embedding(input_file):
    umap_embedding = {}
    with open(input_file, 'r')as f:
        lines = f.readlines()
        for line in lines[1:]:
            line_list = line.split('\t')
            label = line_list[2]
            umap1 = float(line_list[4])
            umap2 = float(line_list[5].replace('\n', ''))
            umap_embedding[label] = [umap1, umap2]
    return umap_embedding

def check_cluster(umap_embedding, color_list):
    for label, embedding in umap_embedding.items():
        if embedding[0]<11:
            print(f'CONTAINS=={label} {color_list[0]}')
        elif embedding[0]>27:
            print(f'CONTAINS=={label} {color_list[1]}')
        else:
            print("out the range")
umap_embedding = read_embedding(input_file='./result/all_site_result/umap/species/Enneapterygius_etheostomus/sequences.tsv')
check_cluster(umap_embedding=umap_embedding, color_list=['#0000ff', '#00ff00'])