# Global Ocean Beaching Pathways of Surface Drifters

This repository analyzes beaching patterns of ocean surface drifters from the [Global Drifter Program (GDP)](https://www.aoml.noaa.gov/global-drifter-program/). We use observational drifter trajectories to understand surface ocean circulation and coastal beaching mechanisms.

---

## Quick Start: Key Figures & Workflows

### 1. Data Description Histograms
**Notebook:** `data_stats_histograms.ipynb`

Explore the dataset characteristics:
- How long it takes drifters to beach
- How far they travel

**Output:** Statistical plots and histograms of raw data

---

### 2. Beaching Probability Heatmaps
**Notebook:** `probability_of_beaching.ipynb`

Generate probability maps showing ocean regions that produce higher-probability-of-beaching drifters:
- Global beaching likelihood by grid cell
- Regional probability distributions

**Output:** Heatmaps of beaching probability across the global ocean

---

### 3. Time to Beaching Heatmaps
**Notebook:** `time_to_beach.ipynb`

Calculate how long drifters take to reach the coast:
- Expected travel time from origin to beaching
- Temporal patterns by region
- Distribution of beaching timescales

**Output:** Heatmaps of time to beach

---

### 4. Watershed Delineation & Density Maps
**Notebook:** `watershed_figures.ipynb`

Define coastal watersheds and quantify drifter beaching:
- Global ocean delineated by coastal watershed boundaries
- Coastline identification and mapping
- Drifter count density by watershed
- Normalized min/max density visualization

**Output:** Watershed maps and density heatmaps

---

### 5. Watershed Surface Areas
**Notebook:** `watershed_surface_areas.ipynb`

Calculate geometric properties:
- Surface area of each identified watershed
- Area-normalized beaching statistics
- Comparison of relative watershed sizes

**Output:** Watershed area statistics

---

### 6. Wind-Dependent Beaching Probability
**Script:** `plotangles.m` (MATLAB)

Analyze beaching probability as a function of wind:
- Wind direction relative to coastline normal vector
- Wind speed effects on beaching likelihood
- Directional beaching patterns

**Output:** Polar plots and directional probability distributions

---

## Data Processing Pipeline

### Data Source
- **Original data access location:** [GDP Hourly Data](https://www.aoml.noaa.gov/phod/gdp/hourly_data.php)
- **Current location:** `ftp://ftp.aoml.noaa.gov/phod/pub/buoydata/hourly_product/`

### Key Processing Steps

**1. Data Acquisition & Cleaning**
- Download hourly drifter trajectories from GDP
- Filter for quality-controlled observations
- Handle missing data and gaps

**2. Drogue Methods**
- Drifter trajectories divided into drogue-on and drogue-off periods
- Separate into drogued and undrogued datasets
- Further analysis is mostly performed on the undrogued dataset as these can act as better proxies for surface advecting plastics and other objects

**3. Beaching Detection**
- Classify beaching vs. non-beaching trajectories by drifter death code
- Extract beaching location and timing

**4. Probability Calculation**
- Grid-based probability estimation
- Kernel density estimation

**5. Time-to-Beach Analysis**
- Calculate duration from specified origin to beaching point
- Account for trajectory gaps and sampling intervals
- Statistical distributions by region

**6. Clustering Methods**
- Identify geographic clusters of beaching events with ensembles of HDBSCAN and agglomerative clustering experiments
- Hyperparameter sweeps and sensitivity analysis

**7. Coastal Angle & Wind Analysis**
- Compute local coastline normal vectors
- Relate drifter beaching probability to wind direction/speed relative to coast

---

## Repository Structure

```
├── README.md                          # This file
├── data_stats_histograms.ipynb        # Dataset overview & statistics
├── probability_of_beaching.ipynb      # Beaching probability maps
├── time_to_beach.ipynb                # Time-to-beaching analysis
├── watershed_figures.ipynb            # Watershed maps & density
├── watershed_surface_areas.ipynb      # Watershed area calculations
├── plotangles.m                       # Wind-dependent beaching analysis
├── data_processing/                   # Data handling scripts
├── drogue_analysis/                   # Drogue-specific analysis
├── windland_analysis/                 # Wind & land analysis
├── watersheds/                        # Watershed data & outputs
├── ML/                                # Some very old machine learning experiments (predict trajectories, beaching coordinates)
├── old/                               # Even older notebooks/scripts/experiments
└── LICENSE                            # MIT License
```

---

## Dependencies

- Python 3.x with Jupyter Notebook
- Scientific Python stack: NumPy, SciPy, Pandas, Matplotlib
- Geospatial tools: Cartopy, Basemap (or similar)
- MATLAB (for `plotangles.m`)

---

## References

- [Global Drifter Program](https://www.aoml.noaa.gov/global-drifter-program/)
- [GDP Hourly Data Portal](https://www.aoml.noaa.gov/phod/gdp/hourly_data.php)
- [FTP Data Repository](ftp://ftp.aoml.noaa.gov/phod/pub/buoydata/hourly_product/)

---

**License:** MIT License
