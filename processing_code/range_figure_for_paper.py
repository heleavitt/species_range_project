import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature

sns.set(style="white")

# ---- Load data
df = pd.read_csv(

    r'C:\Users\hl51981\OneDrive - University of Georgia\Leavitt_Herbert\PFFW\Manuscripts\Global Change\Revision_repository\outputs\species_latitudes.csv'
)
# Auto-rename XYYYY -> "Year YYYY"
df.columns = [f"Year {c[1:]}" if c.startswith("X") and c[1:].isdigit() else c for c in df.columns]

# ---- Config: column names in your CSV
COL_SPECIES = 'species'
COL_MIN = 'Min Latitude'
COL_MAX = 'Max Latitude'
COL_MEAN = 'Mean Latitude'

# ---- Build species→longitude map, ordered by mean latitude (high -> low)
def build_species_longitudes(d, center_lon=-80, span_deg=90):
    means = (
        d[[COL_SPECIES, COL_MEAN]]
        .dropna()
        .drop_duplicates()
        .assign(spec=lambda x: x[COL_SPECIES].astype(str).str.strip().str.lower())
        .sort_values(COL_MEAN, ascending=False)  # highest first
        .reset_index(drop=True)
    )
    if len(means) == 1:
        lons = np.array([center_lon])
    else:
        lons = np.linspace(center_lon - span_deg/2, center_lon + span_deg/2, len(means))
    return dict(zip(means['spec'], lons))

species_to_lon = build_species_longitudes(df, center_lon=-80, span_deg=90)

# ---- Helper: get set of species present in a given year
def species_in_year_set(d, year_label):
    if year_label not in d.columns:
        raise KeyError(f"Missing column: {year_label}")
    vals = d[year_label]
    if vals.dtype == bool:
        mask = vals
    else:
        mask = vals.astype(str).str.strip().str.lower().isin({'1','true','t','yes','y'})
    return set(d.loc[mask, COL_SPECIES].astype(str).str.strip().str.lower())

# ---- Plotter
def plot_species_latitude_ranges_for_year(ax, year_label, d, species_to_lon):
    present = species_in_year_set(d, year_label)

    # Map
    ax.set_extent([-130, -30, -60, 60], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor='tan')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.6)
    ax.add_feature(cfeature.BORDERS, linestyle=':', linewidth=0.5)
    gl = ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.4,
                      xlabel_style={'size': 10},   # longitude label font size
                      ylabel_style={'size': 10})    # latitude label font size)
    gl.right_labels = False
    gl.top_labels = False


    # Ranges (staggered by species-specific longitude)
    for _, row in d.iterrows():
        sp = str(row[COL_SPECIES]).strip().lower()
        if sp not in species_to_lon:
            continue
        lon = species_to_lon[sp]
        lat_min = float(row[COL_MIN])
        lat_max = float(row[COL_MAX])
        lat_mean = float(row[COL_MEAN])

        in_year = sp in present
        color = 'royalblue' if in_year else 'lightgray'
        mean_color = 'red' if in_year else 'lightgray'

        ax.plot([lon, lon], [lat_min, lat_max], color=color, linewidth=2,
                solid_capstyle='butt', transform=ccrs.PlateCarree())
        ax.plot(lon, lat_mean, marker='o', color=mean_color, markersize=2,
                transform=ccrs.PlateCarree())

    # Reference latitude lines
    for lat in [0, 23.5, -23.5, 35, -35]:
        ax.plot([-130, -30], [lat, lat], linestyle='--', linewidth=0.5, color='k',
                transform=ccrs.PlateCarree())

    ax.set_title(f'Species Latitude Range and\nMean Latitude in {year_label}', fontsize=14)

# ---- Figureimport matplotlib.gridspec as gridspec

# ---- Figure: 18 x 24 cm
fig = plt.figure(figsize=(7.09, 9.45))  # 18 cm wide x 24 cm tall
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])  # 2 rows x 2 cols

# Top axis (spans both columns)
ax_top = fig.add_subplot(gs[0, :], projection=ccrs.PlateCarree())

# Bottom row: two side-by-side
ax_bl = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
ax_br = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree())

# ---- Plots
plot_species_latitude_ranges_for_year(ax_top, 'Year 2005', df, species_to_lon)
plot_species_latitude_ranges_for_year(ax_bl, 'Year 2016', df, species_to_lon)
plot_species_latitude_ranges_for_year(ax_br, 'Year 2022', df, species_to_lon)

# Adjust spacing
fig.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.01, wspace=0.15, hspace=0.1)

# Save
plt.savefig(
    r'C:\Users\hl51981\OneDrive - University of Georgia\Leavitt_Herbert\PFFW\Manuscripts\Global Change\Revision_repository\outputs\plots\species_range_fig.png',
    dpi=600
)
plt.show()
