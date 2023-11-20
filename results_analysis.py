import os
import sys
import pandas as pd
import matplotlib.pyplot as plt


def define_arguments(): # Fetch and format arguments

    if len(sys.argv) != 4: # Check that 4 arguments are given
        print("This script requires 4 inputs: path/to/script, path/to/trimmed/reads, path/to/intersected/reads, path/to/AsiSI/cut/sites")
        sys.exit(1) # Exit process and raise a SystemExit exception
    
    trimmed_input_dir = sys.argv[1]
    if not trimmed_input_dir.endswith('/'): # Ensure the trimmed input directory path ends with a slash
        trimmed_input_dir += '/'

    intersected_input_dir = sys.argv[2] # Assign the 3rd argument as the intersected .bed file directory
    if not intersected_input_dir.endswith('/'): # Ensure the intersected input directory path ends with a slash
        intersected_input_dir += '/'

    asiSI_cut_sites = sys.argv[3]

    return trimmed_input_dir, intersected_input_dir, asiSI_cut_sites


def find_total_breaks(intersected_input_dir):

    intersected_columns = ['chromosome_id', 'start', 'end', 'break_num', 'q_score', 'strand'] # Define the dataframe column names

    sample_and_total_breaks_list = [] # Define the temporary master list

    for file in os.listdir(intersected_input_dir):
        if file.endswith('.bed'): # Only capture .bed files
            bed_df = pd.read_csv(f'{intersected_input_dir}{file}', sep='\t', header=None, names=intersected_columns) # Read the tab delimited file and save to a pandas dataframe

            total_breaks = bed_df['break_num'].sum() # Final the total number of breaks for the given sample
            sample_id = file.replace('intersected_trimmed_', '') # Remove the file prefixes
            sample_id = sample_id.replace('.bed', '') # Remove the file suffix

            data_list = [sample_id, total_breaks] # Assign the sample_id and total_breaks to a temporary list
            sample_and_total_breaks_list.append(data_list) # Append temporary list to the master list

    df1 = pd.DataFrame(sample_and_total_breaks_list, columns=['sample_id', 'total_breaks']) # Create pandas dataframe from the master list

    return df1


def find_row_count(trimmed_input_dir):

    trimmed_columns = ['chromosome_id', 'start', 'end', 'id', 'q_score', 'strand'] # Define the dataframe column names

    sample_and_row_count_list = [] # Define the temporary master list

    for trimmed_file in os.listdir(trimmed_input_dir):
        if trimmed_file.endswith('.bed'): # Only capture .bed files
            bed_df = pd.read_csv(f'{trimmed_input_dir}{trimmed_file}', sep='\t', header=None, names=trimmed_columns) # Read the tab delimited file and save to a pandas dataframe
            row_count = len(bed_df) # Count the number of rows in the dataframe
            row_count = row_count/1000 # Divide row count by 1000

            sample_id = trimmed_file.replace('trimmed_', '') # Remove the file prefix
            sample_id = sample_id.replace('.bed', '') # Remove the file suffix

            data_list = [sample_id, row_count] # Assign the sample_id and total_breaks to a temporary list
            sample_and_row_count_list.append(data_list) # Append temporary list to the master list

    df2 = pd.DataFrame(sample_and_row_count_list, columns=['sample_id', 'row_count/1000'])  # Create pandas dataframe from the master list

    return df2


def find_cut_site_percentage(asiSI_cut_sites, intersected_input_dir):

    asiSI_sites_df = pd.read_csv(asiSI_cut_sites, sep='\t', header=None)
    total_asiSI_sites = len(asiSI_sites_df)

    intersected_columns = ['chromosome_id', 'start', 'end', 'break_num', 'q_score', 'strand']

    sample_and_nonzero_breaks_list = []  # Initialize an empty list to store sample_id and non-zero breaks

    for file in os.listdir(intersected_input_dir):
        if file.endswith('.bed'):  # Only capture .bed files
            bed_df = pd.read_csv(f'{intersected_input_dir}{file}', sep='\t', header=None, names=intersected_columns)
            
            non_zero_breaks = (bed_df['break_num'] != 0).sum()  # Count non-zero values in the 'break_num' column

            asiSI_cut_percent = (non_zero_breaks/total_asiSI_sites) * 100 # Find the percentage number of AsiSI cut sites

            sample_id = file.replace('intersected_trimmed_', '').replace('.bed', '')  # Remove file prefixes and suffix
            
            data_list = [sample_id, asiSI_cut_percent]  # Assign sample_id and non-zero breaks to a temporary list
            sample_and_nonzero_breaks_list.append(data_list)  # Append temporary list to the master list

    df3 = pd.DataFrame(sample_and_nonzero_breaks_list, columns=['sample_id', 'asiSI_cut_percent'])

    return df3


def merge_dataframes(df1, df2, df3):

    merged_df = pd.merge(df1, df2, on='sample_id') # Merge the two dataframes on the sample_id
    merged_df = pd.merge(merged_df, df3, on='sample_id')

    merged_df['normalised_break_count'] = merged_df['total_breaks'] / merged_df['row_count/1000'] # Calculate normalised break count
    
    # Preferential sample order
    ascending_order = ['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5', 'Sample6', 'Sample7', 'Sample8', 'Sample9', 'Sample10', 'Sample11', 'Sample12', 'Sample13', 'Sample14', 'Sample15', 'Sample16']

    merged_df = merged_df.set_index('sample_id').loc[ascending_order].reset_index() # Reorder merged_df based on the prefered order

    print(merged_df)

    return merged_df


def plot_data(merged_df):

    plt.figure(figsize=(8, 8)) # Define the plot dimensions
    plt.scatter(merged_df['sample_id'], merged_df['normalised_break_count'], color='black') # Create a scatter plot

    plt.xlabel('Sample') # Set the x axis label
    plt.xticks(rotation=90, fontsize=7) # Rotate and resize x axis text
    plt.ylabel('Normalised Break Count') # Set the y axis label
    #plt.ylim(0,5)

    plt.savefig('scatter_plot.png') # Save the plot


def main(): # Script workflow

    trimmed_input_dir, intersected_input_dir, asi_SI_cut_sites = define_arguments() # Call the define_arguments() function

    df1 = find_total_breaks(intersected_input_dir) # Call the find_total_breaks() function

    df2 = find_row_count(trimmed_input_dir) # Call the find_row_count() function

    df3 = find_cut_site_percentage(asi_SI_cut_sites, intersected_input_dir) # Call the find_cut_site_percentage() function

    merged_df = merge_dataframes(df1, df2, df3) # Call the merge_dataframes() function

    plot_data(merged_df) # Call the plot_data() function


if __name__ == '__main__': # Safety check to only execute the script if the file is run directly

    main() # Call main() function to execute script workflow
