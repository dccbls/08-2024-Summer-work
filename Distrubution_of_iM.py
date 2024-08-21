# %%
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

# Step 1: Load the data from the Excel file
data_filtered_mRNA = pd.read_excel('/Users/yding/Desktop/filter_CS_mRNA_conformation.xlsx')

# Step 2: Filter out rows where any of the UTR or CDS lengths are zero
data_filtered_mRNA = data_filtered_mRNA[(data_filtered_mRNA['5UTR_length'] > 0) & 
                                        (data_filtered_mRNA['CDS_Length'] > 0) & 
                                        (data_filtered_mRNA['3UTR_Length'] > 0)]

# Step 3: Calculate the relative position across the entire mRNA
data_filtered_mRNA['relative_position'] = data_filtered_mRNA.apply(
    lambda row: row['start'] / row['mRNA_Length'], axis=1)

# Step 4: Calculate the total number of unique genes
total_genes = data_filtered_mRNA['ID'].nunique()

# Step 5: Define the binning for the entire mRNA
bins_mRNA = 100

# Step 6: Calculate histogram for the entire mRNA
counts_mRNA, bin_edges = np.histogram(data_filtered_mRNA['relative_position'], bins=bins_mRNA)

# Step 7: Normalize the counts by the total number of genes
counts_mRNA = counts_mRNA / total_genes

# Step 8: Generate x values based on the bin edges
x_values = (bin_edges[:-1] + bin_edges[1:]) / 2 * 100  # Convert to percentage

# Step 9: Annotate regions based on x_values
utr_5_length_percentage = data_filtered_mRNA['5UTR_length'].mean() / data_filtered_mRNA['mRNA_Length'].mean() * 100
utr_3_length_percentage = data_filtered_mRNA['3UTR_Length'].mean() / data_filtered_mRNA['mRNA_Length'].mean() * 100

region_labels = []
for x in x_values:
    if x <= utr_5_length_percentage:
        region_labels.append("5' UTR")
    elif x <= 100 - utr_3_length_percentage:
        region_labels.append('CDS')
    else:
        region_labels.append("3' UTR")

# Adding the region labels to the DataFrame
output_data = pd.DataFrame({
    'x_values': x_values,
    'counts_mRNA': counts_mRNA,
    'Region': region_labels
})

# Save the annotated data to a CSV file
output_file_path_annotated = '/Users/yding/Desktop/i_motif_distribution.csv'
output_data.to_csv(output_file_path_annotated, index=False)

# Step 10: Smooth the curve using spline interpolation
x_smooth = np.linspace(x_values.min(), x_values.max(), 300)
spl = make_interp_spline(x_values, counts_mRNA, k=3)
y_smooth = spl(x_smooth)

# Step 11: Plotting the data
plt.figure(figsize=(12, 6))

plt.axvspan(0, utr_5_length_percentage, color='orange', alpha=0.3, label="5' UTR")
plt.axvspan(utr_5_length_percentage, 100 - utr_3_length_percentage, color='gray', alpha=0.3, label='CDS')
plt.axvspan(100 - utr_3_length_percentage, 100, color='purple', alpha=0.3, label="3' UTR")

# Plot the smoothed data
plt.plot(x_smooth, y_smooth, color='black')

plt.title('Normalized Distribution of i-motif Structures Across Gene Regions')
plt.xlabel('Relative Position (Percentage)')
plt.ylabel('Average Number of i-motif Structures per Gene')
plt.legend(loc='upper right')
plt.grid(True)

# Show the plot
plt.show()

# Output the file path for CSV file download
output_file_path_annotated


