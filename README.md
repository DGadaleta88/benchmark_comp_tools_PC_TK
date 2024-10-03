# Project Name: Benchmarking of Computational Tools for Physicochemical and Toxicokinetic Property Predictions

This repository contains the code and datasets used for benchmarking various computational tools for predicting physicochemical (PC) and toxicokinetic (TK) properties of chemicals. The work presented here evaluates the performance of different freely available tools, aiming to provide useful recommendations for researchers, regulatory authorities, and industries.

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [File description](#file-description)
4. [Data Description](#data-description)
5. [Reference](#reference)
6. [License](#license)

## Overview
This project aims to benchmark several well-known computational tools for the prediction of key chemical properties:
- Physicochemical properties: logP, water solubility, melting point, boiling point, henry's law constant, vapor pressure, pka
- Toxicokinetic properties: plasma protein binding, blood-brain barrier permeability, Caco2 permeability, skin permeation, bioavailability, human intestinal absorption

The key goal is to identify the most reliable tools for high-throughput assessment of these properties, enabling their use in regulatory and industrial settings. The tools that have been benchmarked are:
- OPERA v2.9, VEGA v1.2.3, Episuite v4.1, TEST v5.1.2, ProtoPRED v1.0, ADMETlab 3.0, OCHEM v4.3.156, ChemAxon Calculator Plugin v2.4.1.0, SwissADME, vNN-ADMET, pkCSM, BAYER in-house models

## Installation
To run the code provided in this repository, you need the following dependencies:

- Python 3.8+
- RDKit
- NumPy
- Pandas
- BayesIO
- Matplotlib
- Seaborn
- Math

## File dscription
The files in this repository were used to perform statistical analysis and to generate datasets and figures included int he paper "Comprehensive Benchmarking of Computational Tools for Predicting Toxicokinetic and Physicochemical Properties of Chemicals" (Submitted manuscript).

1. The file benchmark_comp_tools_p1.py was used to geenrate curated datasets used for benchmarking and statistical analysis
2. The file benchmark_comp_tools_p2.knwf was used for generating Additional File 3 and Figures 2, 3, 4 and 5
   
The datasets and the figures generated are available in the data/ folder. These include the curated chemical datasets in CSV format, along with the associated experimental property values.

## Data Description
The datasets used for this benchmarking are collected from multiple literature sources and divided into categories such as drugs, natural products, and industrial chemicals. These datasets have been carefully curated to remove structural outliers, and experimental values have been cross-validated across different datasets.
To ensure data quality, the following steps were applied during curation:

1. Response Outlier Removal: Z-scores were calculated for each dataset, and data points with a Z-score > 3 were considered outliers and removed.
2. Cross-dataset Comparison: Experimental property values of chemicals shared among different datasets were compared. Chemicals with inconsistent values (standard deviation > 0.2 across datasets) were removed.
3. Final Dataset: The curated datasets are available in the data/ folder in the present workflow area accessible from your KNIME workspace (C:\Users\<username>\knime-workspace\<workflow-pah>\Additional File 5\data.

Performance of tools in predicting datasets were evaluated with regard of the models’ external predictivity, emphasizing the performance inside of the models inside the’ applicability domain.

## Reference
Additional details are available from:

Gadaleta D., Serrano Candelas E., Ortega Vallbona R., Colombo E., Garcia de Lomana M., Biava G., Aparicio-Sanchez P., Roncaglioni A., Gozalbes R., Benfenati E. Comprehensive Benchmarking of Computational Tools for Predicting Toxicokinetic and Physicochemical Properties of Chemicals. Submitted manuscript.

## License
This project is licensed under the MIT License. See the LICENSE file for details.
