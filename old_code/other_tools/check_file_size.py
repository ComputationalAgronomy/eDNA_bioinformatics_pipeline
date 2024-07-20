import os

folders = os.listdir('./umap_result')
print(folders)
for folder in folders:
    path = f'./umap_result/{folder}/distance.txt'
    file_size = os.path.getsize(path)
    if file_size/1024 > 10:
        print(folder)