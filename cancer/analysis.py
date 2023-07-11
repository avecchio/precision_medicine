import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster import hierarchy
#hierarchy.set_link_color_palette(['c', 'm', 'y', 'k', 'r'])
rnaseq_df = pd.read_csv("test.rnaseq.txt", sep='\t')

plt.figure(figsize=(10,6))
link = linkage(rnaseq_df, method='complete', metric='euclidean')
#c = fcluster(link, t=0.9, criterion='distance')
#print(c)
#dend = dendrogram(link, p=5, truncate_mode='level')

dendrogram(link,
            orientation='top',
            #labels=['A', 'B', 'C', 'D'],
            distance_sort='descending',
            show_leaf_counts=True)
#print(dend)
plt.title('Dendrogram')
plt.xlabel('Samples')
plt.ylabel('Euclidean distances')
plt.savefig('scipy_dendrogram.png')



def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    #print(model.labels_)
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    d = dendrogram(linkage_matrix, labels=model.labels_, **kwargs)
    label_colors = {'0': 'r', '1': 'g', '2': 'b', '3': 'm'}

    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()
    for lbl in xlbls:
        lbl.set_color(label_colors[lbl.get_text()])

    #print(d)

#from sklearn.cluster import AgglomerativeClustering
#import numpy as np
#X = np.array([[1, 2],  [1, 0],
#               [4, 2], [4, 4], [1, 4], [4, 0]])
#clustering = AgglomerativeClustering().fit(X)
#print(clustering.labels_)

# setting distance_threshold=0 ensures we compute the full tree.
model = AgglomerativeClustering(n_clusters=4, affinity="euclidean", linkage="complete", compute_distances=True) #, distance_threshold=None)

model = model.fit(rnaseq_df.transpose())
#print(model.feature_names_in_)
print(model.labels_)
group_map = {}
samples = rnaseq_df.columns
labels = model.labels_
for index, sample in enumerate(samples):
    group_map[sample] = labels[index]

print(group_map)

plt.figure(figsize=(10,6))
plt.title("Hierarchical Clustering Dendrogram")
# plot the top three levels of the dendrogram
plot_dendrogram(model)
#print(model.children_)
plt.xlabel("Number of points in node (or index of point if no parenthesis).")
plt.savefig('scikitlearn_dendrogram.png')


#print(group_map)

chr = 'chr17'
start = 37844347
end = 37884911
#chr17:37,844,347-37,884,911
group_count = {}
for label in labels:
    group_count[str(label)] = 0

cnvs_df = pd.read_csv("test.CNV.genomicSegment.txt", sep='\t')
for record in cnvs_df.to_dict('records'):
    if ((record['chr'] == chr)) and ((record['start'] <= int(start)) or (record['end'] >= int(end))):
        if (record['CNV'] > 0):
            sample = record['Sample']
            group_id = group_map[sample]
            group_count[str(label)] += 1
print(group_count)

mutations_df = pd.read_csv("test.mutation.txt", sep='\t')

group_gene_count = {}
for label in labels:
    group_gene_count[label] = {}
for record in mutations_df.to_dict('records'):
    sample = record['sample']
    gene = record['gene']

    group_id = group_map[sample]
    if gene not in group_gene_count[label]:
        group_gene_count[label][gene] = 0
    group_gene_count[label][gene] += 1
#print(group_gene_count)
