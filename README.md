# GPB Field Data Analysis Shiny App

**Author:** Udaya Bhanu Angirekula  
**Email:** ubhanu2698@gmail.com  
**GitHub:** [https://github.com/ubhanu2698/GPB_Field_Data_Analysis](https://github.com/ubhanu2698/GPB_Field_Data_Analysis)  
**Shiny App:** [https://gpbtools.shinyapps.io/GPB_feild_data_analysis_v1/](https://gpbtools.shinyapps.io/GPB_feild_data_analysis_v1/)

---

## Overview

The **GPB Field Data Analysis Shiny App** is an interactive web application designed for analyzing multi-environment wheat field trial data.  
It allows researchers and breeders to perform statistical analyses, including ANOVA, BLUP/BLUE estimation, and data visualization, across multiple environments, years, and replications.  

This app was developed using **R** and **Shiny** and is hosted online via **shinyapps.io** for easy access.

---

## Features

- Upload field trial datasets in **CSV** or **Excel** format  
- Perform **ANOVA** for multi-environment trials  
- Calculate **BLUPs and BLUEs** for genotypes  
- Visualize genotype performance across environments  
- Supports multiple years, environments, replications, and blocks  
- Sample dataset included for testing and demonstration  

---

## Installation

1. Clone this repository:

```bash
git clone https://github.com/ubhanu2698/GPB_Field_Data_Analysis.git
```

2. Install required R packages:

```R
install.packages(c("shiny", "readxl", "dplyr", "ggplot2", "PBTools"))
```

3. Run the app locally:

```R
library(shiny)
runApp("path_to_app_folder")
```

---

## Usage

1. Upload your dataset using the file input panel.  
2. Select the analysis type (ANOVA, BLUP/BLUE).  
3. Visualize and download results as needed.  

---

## Sample Data

A small sample dataset (`Sample_Data.xlsx`) is included in this repository for testing the app.

---

## Citation

If you use this app in your research, please cite:  

> Udaya Bhanu Angirekula (2025). GPB Field Data Analysis Shiny App. GitHub repository, [https://github.com/ubhanu2698/GPB_Field_Data_Analysis](https://github.com/ubhanu2698/GPB_Field_Data_Analysis). DOI to be added once archived on Zenodo.

---

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.
