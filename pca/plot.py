import sys

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.markers import MarkerStyle

sns.set_style('ticks')

experiment_dir = sys.argv[1]

coeff_data = pd.read_table(f'/data/projects/varlock/pca/{experiment_dir}/merged_samples_snp.eigenval',
                           sep='\t',
                           index_col=None,
                           header=None)
pc_data = pd.read_table(f'/data/projects/varlock/pca/{experiment_dir}/merged_samples_snp.eigenvec',
                        sep='\t',
                        index_col='IID')

print(coeff_data[0].round(2).tolist())
del pc_data['FID']
all_data = pc_data[['PC1', 'PC2']]
original_data = all_data[~all_data.index.str.endswith('masked')]
masked_data = all_data[all_data.index.str.endswith('masked')]

# print(original_data)
# print(masked_data)
# print(pc_data)
# exit(0)

fig, ax = plt.subplots()

labels = ['original', 'masked']
colors = ['navy', 'darkorange']

for sample in original_data.index:
    x1 = original_data.loc[sample][0]
    y1 = original_data.loc[sample][1]
    x2 = all_data.loc[sample + '_masked'][0]
    y2 = all_data.loc[sample + '_masked'][1]
    ax.arrow(x1, y1, x2 - x1, y2 - y1,
             alpha=0.1, width=0.0001, length_includes_head=True, shape='full', overhang=0, fc='k', ec='k')

ax.scatter(
    original_data.iloc[:, 0],
    original_data.iloc[:, 1],
    edgecolors='black',
    alpha=0.5,
    s=20,
    linewidths=1,
    label=labels[0],
    marker=MarkerStyle(marker='o', fillstyle='none'),
    facecolors='none'
)
ax.scatter(
    masked_data.iloc[:, 0],
    masked_data.iloc[:, 1],
    edgecolors='black',
    alpha=0.5,
    s=20,
    linewidths=1,
    label=labels[1],
    marker=MarkerStyle(marker='v', fillstyle='none'),
    facecolors='none'
)

# percentage variance explained
pve_pc1 = coeff_data[0][0] / sum(coeff_data[0]) * 100
pve_pc2 = coeff_data[0][1] / sum(coeff_data[0]) * 100

plt.xlabel(f'PC1 ({pve_pc1:.2f} %)')
plt.ylabel(f'PC2 ({pve_pc2:.2f} %)')

plt.tick_params(axis='both', which='both', bottom=False, left=False, labelbottom=False, labelleft=False)

# plt.xlim((-0.1, 0.1))
# plt.ylim((-0.1, 0.1))

# for sample in all_data.index:
#     x = all_data.loc[sample][0]
#     y = all_data.loc[sample][1]
#     text = ax.annotate(sample, xy=(x, y), xytext=(2, 2), textcoords='offset points')
#     text.set_alpha(.5)

plt.legend(loc='best', shadow=False, scatterpoints=1)
# plt.title('PCA of genotypes')
plt.savefig(f'/data/projects/varlock/pca/{experiment_dir}/result.png', dpi=300, pad_inches=0)

plt.show()
