# Monthly heat-stress indicators from ERA5 and CMIP6

This repository provides a fully reproducible and FAIR-oriented pipeline to compute monthly heat-stress indicators from ERA5 reanalysis and CMIP6 climate projections (NASA NEX-GDDP-CMIP6).

The workflow is designed to support climate impact analyses focused on:

human habitability constraints,

outdoor work capacity limitations,

and inter-comparison of commonly used heat-stress indices.

All computations are performed at daily resolution, then aggregated to monthly metrics, which constitute the only final outputs.

Indicators computed

The pipeline computes the following monthly indicators:

Core impact indicators

Wet-bulb temperature (Tw)
Stull (2011), computed from daily mean temperature and relative humidity
→ monthly mean
→ monthly number of days with Tw ≥ 35 °C (habitability threshold)

Wet-Bulb Globe Temperature (WBGT, proxy – shade/indoor)
WBGT = 0.7 × Tw + 0.3 × Tmean
→ monthly mean
→ monthly number of days with WBGT ≥ 30 °C and WBGT ≥ 32 °C
(work-capacity thresholds)

Complementary indicators (comparison & robustness)

Heat Index (HI) (NOAA/Rothfusz formulation, using daily maximum temperature)
→ monthly mean only

Humidex (Canadian formulation)
→ monthly mean only

No threshold-based impact metrics are derived from HI or Humidex; they are included to facilitate comparison with previous studies and to assess inter-indicator consistency.

Temporal coverage
ERA5 (reference)

1991–2020
→ 360 monthly values per variable

CMIP6 (NEX-GDDP)

Historical: 1991–2014

SSP1-2.6: 2015–2100

SSP2-4.5 & SSP5-8.5: 2020–2100

Derived analysis windows:

Pseudo-historical reference: 1991–2020
(historical 1991–2014 + SSP1-2.6 2015–2020)

Mid-century: 2041–2070

End-century: 2071–2100

Each analysis window contains 360 monthly values (30 years × 12 months).

Input data structure

This pipeline intentionally preserves the original directory and file naming conventions of the NASA NEX-GDDP-CMIP6 archive.

The expected input structure is identical to the official distribution:

NEX-GDDP/
 └── MODEL/
     └── SCENARIO/
         └── VARIABLE/
             └── VARIABLE_day_MODEL_SCENARIO_*_YEAR_v2.0.nc


where * corresponds to the model-specific ensemble member and grid identifier
(e.g. r1i1p1f1_gr1, which may vary across models).

This design choice ensures:

full compatibility with the official NASA archive,

minimal data duplication,

transparent data provenance,

straightforward reuse of existing local mirrors.

ERA5 inputs follow the standard ECMWF NetCDF daily format.

No original ERA5 or CMIP6 data are modified, redistributed, or duplicated.

Outputs

Final NetCDF outputs contain only aggregated monthly variables:

monthly means of all indicators,

monthly counts of threshold exceedances (Tw, WBGT).

All intermediate daily variables are discarded.

Outputs are organised into:

scenario-level monthly files,

analysis-ready “packs” corresponding to each 30-year study window.

These files are intended to be used directly for figure generation and spatial analysis (e.g. in R).

Workflow orchestration

The pipeline is orchestrated using Snakemake, which:

manages dependencies automatically,

parallelises computations across models and scenarios,

allows safe restart and partial recomputation.

Typical execution:

snakemake -j 4

Dependencies

CDO (Climate Data Operators)

Snakemake

Unix shell environment (tested under WSL/Linux)

Data availability

Input data must be obtained separately from:

ERA5 (ECMWF)

NASA NEX-GDDP-CMIP6

This repository contains only code and configuration files, in accordance with data licensing and redistribution policies.
