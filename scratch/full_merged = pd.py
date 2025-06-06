full_merged = pd.merge(
    species_codex,
    unique_taxa_df,
    left_on='full_name ',
    right_on='Taxon',
    how='outer'
)

full_merged['full_name '] = full_merged['full_name '].fillna(full_merged['Taxon'])
output_path = "full_merged_updated.csv"


full_merged.to_csv(output_path, index=False)
# manual editng happened here

full_merge_edited = pd.read_csv("full_merged_updated_edited.csv")

species_codes= pd.merge(study_species, full_merge_edited, left_on = 'valid_name', right_on = 'full_name ',     how='outer')
species_codes.to_csv("final_species_codex.csv")

# Get unique values in the 'taxon' column
unique_taxa = merged_df['taxon'].unique()