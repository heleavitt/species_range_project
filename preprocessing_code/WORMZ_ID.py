import pandas as pd
import requests
import os
os.chdir('C:\\Users\\hl51981\\OneDrive - University of Georgia\\Leavitt_Herbert\\PFFW\\Manuscripts\\Global Change\\Revision_repository')

def get_worms_id_by_scientific_name(scientific_name):
    scientific_name = str(scientific_name)  # Ensure it's a string
    url = "https://www.marinespecies.org/rest/AphiaRecordsByName/" + scientific_name
    params = {
        'like': 'false',
        'offset': '1'
    }
    response = requests.get(url, params=params, headers={'accept': 'application/json'})
    
    if response.status_code == 200:
        data = response.json()
        if data:
            first_result = data[0]
            worms_id = first_result.get('AphiaID')
            true_name = first_result.get('valid_name', first_result.get('scientificname'))
            rank = first_result.get('rank', 'Unknown')
            return worms_id, true_name, rank
        else:
            return None, "No results found", "Unknown"
    else:
        print(f"{scientific_name} API request failed with status code {response.status_code}: {response.text}")
        return None, "API request failed", "Unknown"

def process_table(input_file_path, output_file_path):
    # Load the table of scientific names
    df = pd.read_csv(input_file_path)
    
   # Add new columns if they don't exist
    if 'AphiaID' not in df.columns:
        df['AphiaID'] = None
    if 'true_name' not in df.columns:
        df['true_name'] = None
    if 'taxonomic_rank' not in df.columns:
        df['taxonomic_rank'] = None

    for index, row in df.iterrows():
        if pd.isnull(row['AphiaID']) or pd.isnull(row['true_name']) or pd.isnull(row['taxonomic_rank']):

            scientific_name = row['valid_name']
            worms_id, valid_name, rank = get_worms_id_by_scientific_name(scientific_name)
            df.at[index, 'AphiaID'] = worms_id
            df.at[index, 'true_name'] = valid_name
            df.at[index, 'taxonomic_rank'] = rank

    df.to_csv(output_file_path, index=False)
    print(f"Results saved to {output_file_path}")

# Specify your input and output file paths
input_file_path = 'final_species_codex.csv'

output_file_path = 'final_aphia_codex.csv'

# Process the table
process_table(input_file_path, output_file_path)


