#For RLK

#This code is used to plot the physical clustering of genes based on their location on the chromosome. The input data is a BED file that #contains information about the genomic coordinates #and names of each gene. The code calculates the midpoint of each gene and then #identifies genes that are physically clustered based on a user-defined threshold. The output is a scatter plot #where each point #represents a physically clustered group of genes. The x-axis shows the genomic position of the midpoint of each gene and the y-axis is an #arbitrary scale. The code also #adds the names of the genes to the plot using annotations.



import pandas as pd
import matplotlib.pyplot as plt

# Load the BED file
bed_df = pd.read_csv('RLK.bed', sep='\t', header=0, names=['chrom', 'start', 'end', 'name'])

# Convert start and end columns to integers
bed_df['start'] = bed_df['start'].astype(int)
bed_df['end'] = bed_df['end'].astype(int)

# Calculate the midpoint of each gene
bed_df['midpoint'] = (bed_df['start'] + bed_df['end']) / 2

# Set a threshold for physical clustering (e.g., 10 kilobases)
threshold = 10000

# Create a boolean mask for genes that are physically clustered
clustered_mask = bed_df.groupby('chrom')['midpoint'].diff().abs() <= threshold

# Plot the clustered genes
fig, ax = plt.subplots(figsize=(10,5))
ax.scatter(bed_df[clustered_mask]['midpoint'], [1]*sum(clustered_mask))
ax.set_xlabel('Genomic position')
ax.set_yticks([])

# Add gene names to the plot
for idx, row in bed_df[clustered_mask].iterrows():
    ax.annotate(row['name'], xy=(row['midpoint'], 1), xytext=(row['midpoint'], 1.01), ha='center', fontsize=8, rotation=90)

plt.show()



import pandas as pd
import matplotlib.pyplot as plt

# Load the BED file
bed_df = pd.read_csv('genes.bed', sep='\t', header=0, names=['chrom', 'start', 'end', 'name'])

# Convert start and end columns to integers
bed_df['start'] = bed_df['start'].astype(int)
bed_df['end'] = bed_df['end'].astype(int)

# Set a threshold for physical clustering (e.g., 10 kilobases)
threshold = 10000

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

    # Plot the clustered genes for the current chromosome
    ax = axes[i]
    ax.scatter(chrom_df[clustered_mask]['midpoint'], [1] * sum(clustered_mask), s=10, alpha=0.8)
    ax.set_yticks([])
    ax.set_title('Chr{}'.format(chrom), fontsize=12, y=-0.1, x=-0.05)

    # Add gene names to the plot for the current chromosome
    for idx, row in chrom_df[clustered_mask].iterrows():
        ax.annotate(row['name'], xy=(row['midpoint'], 1), xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, rotation=90)
        clustered_genes.append(row['name'])

# Set the x-axis label for the bottom-most subplot
axes[-1].set_xlabel('Genomic position', fontsize=12)

plt.show()

# Export the list of clustered genes in a CSV file
clustered_genes_df = pd.DataFrame(clustered_genes, columns=['Clustered Genes'])
print(clustered_genes_df)
clustered_genes_df.to_csv('clustered_genes.csv', index=False)
