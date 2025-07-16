# Historical Observations Data Platform 
Victoria Ford, Héctor Inda Diaz, Vanessa Machuca, Nicole Keeney, Neil Schroeder, Brenna Norris, Ella Belfer, Beth McClenny, Grace DiCecco<br>
Code associated with PIR-19-006

The Historical Observations Data Platform is a cloud-based, historical weather observations data platform to enable California's energy sector access to high-quality, open climate and weather data. The Platform responds to community partner needs in understanding weather and cliamte information including the severity, duration, frequency, and rate of change over time of extreme weather events, as well as supporting projections downscaling efforts. We implement stringent, customized Quality Assurance/Quality Control (QA/QC) procedures in line with international convention, and updates relevant to the energy sector are accurately captures (such as temperature and precipitation extremes, winds, and solar radiation).

This repository contains the code (via Python scripts and Jupyter Notebooks) associated with the full processing pipeline for data ingestion into the Historical Data Platform.

## 🗂️ Repository Structure

```text
historical-obs-platform/
├── data/
├── environment/               # Files for building the computational environment, including a README with further instructions
├── notebooks/                 # Jupyter notebooks for data visualization and analysis 
├── scripts/                   # Data processing code for all steps of the QAQC process 
│   ├── 1_pull_data/           # Scripts for retrieving/scrape network station data from their respective sources 
│   ├── 2_clean_data/          # Scripts for cleaning individual networks to a consistent standard
│   ├── 3_qaqc_data/           # Scripts to QA/QC stations 
│   ├── 4_merge_data/          # Scripts to close out processing, and standardize to hourly timesteps. Data at conclusion have been fully processed.
│   ├── pcluster/              # Code and shell scripts for running QAQC and merge scripts in an AWS pcluster environment 
│   └── tests/                 # Scripts for testing finalized data products         
└── README.md
```

## License

This project is licensed under the GNU GPLv3 - see the [LICENSE](LICENSE) file for details.

## Support

- 📧 **Email**: [ADD EMAIL HERE](email)
- 🐛 **Issues**: [GitHub Issues](LINK ISSUES HERE)
- 💬 **Discussions**: [GitHub Discussions](LINK DISCUSSIONS HERE)

---

## Contributors

[![Contributors](https://contrib.rocks/image?repo=Eagle-Rock-Analytics/historical-obs-platform)](https://github.com/Eagle-Rock-Analytics/historical-obs-platform/graphs/contributors)
