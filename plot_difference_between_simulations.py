# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 10:52:12 2024

@author: virgi
"""
import numpy as np
import matplotlib.pyplot as plt
import os


def plot_difference_histogram(file1, file2, output_folder, output_name):
    
    data1 = np.loadtxt(file1, delimiter='\t', skiprows=1)  # Skip header row
    data2 = np.loadtxt(file2, delimiter='\t', skiprows=1)
    
    print (data1, data2)
    
    if not np.array_equal(data1[:, 0], data2[:, 0]):
        raise ValueError("The x-values in the two files do not match. Ensure both files have the same x-axis bins.")
    
    # Calculate the difference in y-values
    x_values = data1[:, 0]
    y_diff = data1[:, 1] - data2[:, 1]
    
    # Plot histogram as a bar plot
    plt.bar(x_values, y_diff, width=0.8, align='center', color='blue', edgecolor='black')

    
    # Add labels, title, and legend
    plt.xlabel("x")
    plt.ylabel("Difference in Counts")
    plt.title(f"Difference Histogram: {os.path.basename(file1)} - {os.path.basename(file2)}")
   
    
    # Save plot
    plot_filename = os.path.join(output_folder, output_name + '.png')
    plt.savefig(plot_filename)
    plt.show()
    plt.clf()
    plt.close()

    # # Save the x and y_diff data to a new text file
    # diff_data_filename = os.path.join(output_folder, output_name + '.txt')
    # with open(diff_data_filename, 'w') as f:
    #     f.write("x\ty_diff\n")  # Header
    #     for x_val, y_val in zip(x_values, y_diff):
    #         f.write(f"{x_val}\t{y_val}\n")


# Define subfolder paths for the input files and the output

subfolder1 ="C:\\Users\\virgi\\OneDrive\\Documenti\\GitHub\\Modelling-the-phase-separation-of-transcription-factors-on-DNA\\MC\\alfa_0.2"

subfolder2 = "C:\\Users\\virgi\\OneDrive\\Documenti\\GitHub\\Modelling-the-phase-separation-of-transcription-factors-on-DNA\\MC\\alfa_0.2"

output_folder = "C:\\Users\\virgi\\OneDrive\\Documenti\\GitHub\\Modelling-the-phase-separation-of-transcription-factors-on-DNA\\MC\\Simulations_protein_A\\alfa_0.2\\nA_600_n_3000_cluster_histo_Eaa_0_Ead_0_E_ab_0_E_ba_0.txt"



# File names in each subfolder
file1_name = "nA_600_n_3000_cluster_histo_Eaa_0_Ead_0_E_ab_0_E_ba_0.txt"
file2_name = "nA_600_n_3000_cluster_histo_Eaa_0_Ead_1_E_ab_0_E_ba_0.txt"


file1_path = os.path.join(subfolder1, file1_name)
file2_path = os.path.join(subfolder2, file2_name)


output_name = "difference_histogram_prova"


plot_difference_histogram(file1_path, file2_path, output_folder, output_name)




