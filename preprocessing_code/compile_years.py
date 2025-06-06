# Required libraries
import pandas as pd
from pandasgui import show

fish_2006 = pd.read_csv("raw_data/raw_2006_fish.csv")
# Ensure 'Date' is a datetime column
site_data_2006 = pd.read_csv("raw_data/2006_site_data.csv", encoding="latin1")
invert_2006 = pd.read_csv("raw_data/raw_2006_invert.csv")
peneid_2006 = pd.read_csv("raw_data/raw_2006_peneid.csv")
raw_2016 = pd.read_csv("raw_data/raw_2016.csv")
raw_2022 = pd.read_csv("raw_data/raw_2022.csv")
raw_2023 = pd.read_csv("raw_data/raw_2023.csv")
#species_list = pd.read_csv("raw_data/species_list.csv")
species_codex = pd.read_csv("raw_data/species_codex.csv")
study_species = pd.read_csv("species_presence_2006_2016_2022.csv")
# Concatenate the three dataframes
merged_df = pd.concat([fish_2006, invert_2006, peneid_2006], ignore_index=True)
unique_taxa = merged_df['Taxon'].unique()
unique_taxa_df = pd.DataFrame(unique_taxa, columns=['Taxon'])

species_list = pd.read_csv("final_aphia_codex_edited.csv")

# Clean column names
species_list.columns = species_list.columns.str.strip()
species_codex.columns = species_codex.columns.str.strip()

# Process 2022 data
data_2022 = raw_2022.copy()
data_2022["Year"] = 2022
data_2022 = data_2022.melt(id_vars=["site_date_key", "Year"], var_name="species_code", value_name="Count")
data_2022 = data_2022.rename(columns={"site_date_key": "SampleID"})
data_2022 = data_2022.merge(species_list[["species_code", "valid_name"]], on="species_code", how="left")
data_2022 = data_2022.rename(columns={"valid_name": "Taxon"})
data_2022 = data_2022[["Year", "SampleID", "Taxon", "Count"]]

# Process 2023 data using the same structure
data_2023 = raw_2023.copy()
data_2023["Year"] = 2023
data_2023 = data_2023.melt(id_vars=["site_date_key", "Year"], var_name="species_code", value_name="Count")
data_2023 = data_2023.rename(columns={"site_date_key": "SampleID"})
data_2023 = data_2023.merge(species_list[["species_code", "valid_name"]], on="species_code", how="left")
data_2023 = data_2023.rename(columns={"valid_name": "Taxon"})
data_2023 = data_2023[["Year", "SampleID", "Taxon", "Count"]]

# Combine both years
# Get the union of all columns
all_columns = set(data_2022.columns).union(set(data_2023.columns))

# Add any missing columns to each dataframe and fill with 0
for col in all_columns:
	if col not in data_2022.columns:
		data_2022[col] = 0
	if col not in data_2023.columns:
		data_2023[col] = 0


# Reorder columns to be consistent
data_2022 = data_2022[sorted(all_columns)]
data_2023 = data_2023[sorted(all_columns)]

# Concatenate the dataframes
combined_data = pd.concat([data_2022, data_2023], ignore_index=True)

# Process 2016 data and map abbreviated names to full names using species_codex
doerr_to_fullname = species_list.set_index("doerr_name")["valid_name"].to_dict()
data_2016 = raw_2016.copy()
data_2016["Year"] = 2016
data_2016 = data_2016.drop(columns=["date", "bay"], errors="ignore")
data_2016 = data_2016.melt(id_vars=["site_code", "Year"], var_name="Taxon", value_name="Count")
data_2016 = data_2016.rename(columns={"site_code": "SampleID"})
data_2016["Taxon"] = data_2016["Taxon"].replace(doerr_to_fullname)
# Standardize Taxon names for matching
data_2016["Taxon"] = data_2016["Taxon"].str.strip()

# Identify all Palaemonetes entries (species or genus level)
is_palaemonetes = data_2016["Taxon"].str.startswith("Palaemonetes")

# Replace all matching taxa with the genus label
data_2016.loc[is_palaemonetes, "Taxon"] = "Palaemonetes"

# Group and re-sum in case multiple entries now share the same SampleID and Taxon
data_2016 = (
	data_2016.groupby(["Year", "SampleID", "Taxon"], as_index=False)
	.agg({"Count": "sum"})
)

# Ensure 'Date' is a datetime column
site_data_2006["Date"] = pd.to_datetime(site_data_2006["Date"], errors="coerce")

# Filter for GeneralHabitat == 'Marsh' and month == October
marsh_fall_sites = site_data_2006[
	(site_data_2006["GeneralHabitat"].str.strip().str.lower() == "marsh") &
	(site_data_2006["Date"].dt.month == 10)
]

# Extract list of SampleNumbers to use for subsetting data_2006
valid_sample_ids = marsh_fall_sites["SampleNumber"].unique()
# Standardize 2006 data

fish_2006["Year"] = 2006
peneid_2006["Year"] = 2006
invert_2006["Year"] = 2006
invert_2006["Count"] = 1  # Assume presence


minello_to_fullname = species_list.set_index("2005_name")["valid_name"].to_dict()

data_2006 = pd.concat([
	fish_2006[["Year", "SampleNumber", "Taxon", "Count"]],
	peneid_2006[["Year", "SampleNumber", "Taxon", "Count"]],
	invert_2006[["Year", "SampleNumber", "Taxon", "Count"]]
], ignore_index=True)
data_2006 = data_2006.rename(columns={"SampleNumber": "SampleID"})
data_2006["Taxon"] = data_2006["Taxon"].replace(minello_to_fullname)
data_2006 = data_2006[data_2006["SampleID"].isin(valid_sample_ids)]
# Combine all datasets
all_years = pd.concat([data_2006, data_2016, data_2022], ignore_index=True)
all_years = all_years.dropna(subset=["Taxon"])
all_years = all_years[all_years["Count"] > 0]

# === Abundance summary table ===
abundance_summary = (
	all_years.groupby(["Year", "Taxon"], as_index=False)
	.agg(Total_Count=("Count", "sum"))
	.pivot(index="Taxon", columns="Year", values="Total_Count")
	.fillna(0).astype(int).reset_index()
)

# === Sampling effort by year ===
sites_sampled_per_year = pd.DataFrame({
	"Year": [2006, 2016, 2022],
	"Unique_Sites_Sampled": [
		data_2006["SampleID"].nunique(),
		data_2016["SampleID"].nunique(),
		data_2022["SampleID"].nunique()
	]
})

# === Presence summary: number of sites each taxon was detected in ===
presence_summary = (
	all_years.groupby(["Year", "Taxon"])["SampleID"]
	.nunique()
	.reset_index(name="Sites_Present")
)





# === Merge with total number of sites per year ===
presence_summary = presence_summary.merge(
	sites_sampled_per_year, on="Year", how="left"
)

# === Calculate proportion of sites visited with presence ===
presence_summary["Proportion_of_Sites"] = (
	presence_summary["Sites_Present"] / presence_summary["Unique_Sites_Sampled"]
).round(3)
presence_summary.to_csv("presence_summary.csv")

# === Pivot for viewing proportions by taxon and year ===
presence_pivot = (
	presence_summary.pivot(index="Taxon", columns="Year", values="Proportion_of_Sites")
	.fillna(0)
	.reset_index()
)

presence_pivot.to_csv("presence_pivot.csv")


# === Optional: Sampling effort by season ===
def infer_season(sample_id):
	if pd.isna(sample_id) or not isinstance(sample_id, str):
		return "Unknown"
	if any(s in sample_id.lower() for s in ["spr", "apr", "mar"]):
		return "Spring"
	if any(s in sample_id.lower() for s in ["sum", "jun", "jul", "aug"]):
		return "Summer"
	if any(s in sample_id.lower() for s in ["fall", "sep", "oct", "nov"]):
		return "Fall"
	if any(s in sample_id.lower() for s in ["win", "jan", "feb", "dec"]):
		return "Winter"
	return "Unknown"

def extract_season_fallback(sample_id):
	if not isinstance(sample_id, str):
		return "Unknown"
	if any(x in sample_id for x in ["01", "02", "12"]):
		return "Winter"
	if any(x in sample_id for x in ["03", "04", "05"]):
		return "Spring"
	if any(x in sample_id for x in ["06", "07", "08"]):
		return "Summer"
	if any(x in sample_id for x in ["09", "10", "11"]):
		return "Fall"
	return "Unknown"

season_effort = all_years.copy()
season_effort["Season"] = season_effort["SampleID"].apply(infer_season)
season_effort.loc[season_effort["Season"] == "Unknown", "Season"] = (
	season_effort.loc[season_effort["Season"] == "Unknown", "SampleID"]
	.apply(extract_season_fallback)
)

sampling_by_season = (
	season_effort.groupby(["Year", "Season"])["SampleID"]
	.nunique().reset_index(name="Num_Samples")
	.sort_values(["Year", "Season"])
)

# These final tables are:
# - abundance_summary
# - sites_sampled_per_year
# - sampling_by_season
