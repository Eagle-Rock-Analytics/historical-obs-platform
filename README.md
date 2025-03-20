# Historical Observations Data Platform 
Victoria Ford, Héctor Inda Diaz, Beth McClenny, Vanessa Machuca, Nicole Keeney, Brenna Norris, Ella Belfer, Grace DiCecco<br>
Code associated with PIR-19-006

The Historical Observations Data Platform is a cloud-based, historical weather observations data platform to enable California's energy sector access to high-quality, open climate and weather data. The Platform responds to community partner needs in understanding weather and cliamte information including the severity, duration, frequency, and rate of change over time of extreme weather events, as well as supporting projections downscaling efforts. We implement stringent, customized Quality Assurance/Quality Control (QA/QC) procedures in line with international convention, and updates relevant to the energy sector are accurately captures (such as temperature and precipitation extremes, winds, and solar radiation).

This repository contains the code (via Python scripts and Jupyter Notebooks) associated with the full processing pipeline for data ingestion into the Historical Data Platform.
* **1_pull_data**: Scripts to retrieve/scrape wx observation network stations.
* **2_clean_data**: Scripts to clean individual networks to a consistent standard, including metadata, missing values, and data encoding.
* **3_qaqc_data**: Scripts to QA/QC stations, including visualization of tests.
* **4_merge_data**: Scripts to close out processing, and standardize to hourly timesteps. Data at conclusion have been fully processed.

**History**<br>
**v1.0** -- *In progress*
