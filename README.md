# CRISPRi_Screening_AceticAcid

Phenotypic screening of a CRISPRi strain library [Smith et al., 2017](https://doi.org/10.15252/msb.20167233) for acetic acid tolerance

## PURPOSE 
We screened a CRISPR interference library consisting of >9000 Saccharomyces cerevisiae strains where >98% of all essential and respiratory growth-essential genes were targeted with multiple gRNAs. This repository is created to offer any user the opportunity to reproduce the entire analysis starting from the raw data.

## CONTENT 
The codes described in this project shows how to **(a) prepare, (b) analyze, and (c) visualize** the high-throughput phenomics data generated in this study. 

In this project we used two high-thoughput phenomics platform 
* Scan-o-matic [Zackrisson et al., 2016](https://doi.org/10.1534/g3.116.032342)
* Bioscreen microcultivation followed by PRECOG analysis [Fernandez-Ricaud et al, 2016](https://doi.org/10.1186/s12859-016-1134-2) 

All raw data files generated in this study are available in the folder **RAW_DATA**. It has the following sub folders

* **SOM_AA_TITRA**  : Scan-o-matic Acetic acid titration raw data (Plate 7 & 8)
* **SOM_ATC_TITRA** : Scan-o-matic ATc (anhydrous tetracyclin) titration raw data (Plate 7 & 8)
* **SOM_SCR_R001**  : Scan-o-matic CRISPRi library screening raw data, Round 1
* **SOM_SCR_R002**  : Scan-o-matic CRISPRi library screening raw data, Round 2
* **BS_VAL_SCR**    : Bioscreen Validation experiment with selected strains from CRISPRi library raw data.   
              + **STRAIN_MAP_VAL_EXP**: Strains layout in bioscreen plate (well no > strain name).   
              + **BS_PRECOG_OUTPUT**  : PRECOG output (Processed and calibrated growth curve, First derivative, raw phenotypes)

The repository also provides some compiled data for the ease of analysis in the **COMPILED_DATA** folder. 
Those files will be used while running the script in the CRISPRi_Screening_AceticAcid.Rmd file. 

All files in the above folders, the data, and the scripts used to process the data are explained STEP BY STEP in the CRISPRi_Screening_AceticAcid.Rmd file. User can also use the output CRISPRi_Screening_AceticAcid.html file to view the flow and the content of the analysis. 

## CONTRIBUTION
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## LICENSE
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)

