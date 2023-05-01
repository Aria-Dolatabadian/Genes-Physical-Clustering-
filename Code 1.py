#For total genes 

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the BED file
bed_df = pd.read_csv('total genes.bed', sep='\t', header=0, names=['chrom', 'start', 'end', 'name', 'type'])

# Convert start and end columns to integers
bed_df['start'] = bed_df['start'].astype(int)
bed_df['end'] = bed_df['end'].astype(int)

# Set a threshold for physical clustering (e.g., 10 kilobases)
threshold = 10000

# Find the unique gene types in the BED file
gene_types = np.unique(bed_df['type'])

# Assign a different color to each gene type
colors = plt.cm.get_cmap('tab10')(np.arange(len(gene_types)))

# Loop over chromosomes and create subplots
chromosomes = bed_df['chrom'].unique()
num_chromosomes = len(chromosomes)
clustered_genes = []

fig, axes = plt.subplots(nrows=num_chromosomes, figsize=(10, 5 * num_chromosomes), sharex=True)
fig.subplots_adjust(hspace=0.3)

for i, chrom in enumerate(chromosomes):
    # Subset the DataFrame for the current chromosome
    chrom_df = bed_df.loc[bed_df['chrom'] == chrom].copy()

    # Calculate the midpoint of each gene for the current chromosome
    chrom_df.loc[:, 'midpoint'] = (chrom_df['start'] + chrom_df['end']) / 2

    # Create a boolean mask for genes that are physically clustered for the current chromosome
    clustered_mask = chrom_df.groupby('chrom')['midpoint'].diff().abs() <= threshold

    # Plot the clustered genes for the current chromosome, with different colors for each gene type
    ax = axes[i]
    for j, gene_type in enumerate(gene_types):
        gene_type_mask = clustered_mask & (chrom_df['type'] == gene_type)
        ax.scatter(chrom_df[gene_type_mask]['midpoint'], [1] * sum(gene_type_mask), s=10, alpha=0.8, color=colors[j])
    ax.set_yticks([])
    ax.set_title('Chr{}'.format(chrom), fontsize=12, y=-0.1, x=-0.05)

    # Add gene names to the plot for the current chromosome
    for idx, row in chrom_df[clustered_mask].iterrows():
        ax.annotate(row['name'], xy=(row['midpoint'], 1), xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, rotation=90)
        clustered_genes.append(row['name'])

# Set the x-axis label for the bottom-most subplot
axes[-1].set_xlabel('Genomic position', fontsize=12)

# Create a legend showing the gene types and their corresponding colors
handles = [plt.plot([],[], marker="o", ls="", color=color, label=gene_type)[0] for color, gene_type in zip(colors, gene_types)]
plt.legend(handles, gene_types, title='Gene Types', bbox_to_anchor=(1.135, 10), loc='center right', fontsize=10)

plt.show()

# Export the list of clustered genes in a CSV file
clustered_genes_df = pd.DataFrame(clustered_genes, columns=['Clustered Genes'])
print(clustered_genes_df)
clustered_genes_df.to_csv('clustered_genes_total_genes.csv', index=False)
