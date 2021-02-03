CRISPRi\_Phenomics\_Analysis
================
Vaskar Mukherjee
2/3/2021

# IMPORT SCAN-O-MATIC RAW DATA

The phenotypic data generated in scan-o-matic screening in .csv format.
We extract both the absolute and the normalized phenotypes.

The CRISPRi strains in the library were arrayed in 24 plates in 384
format. Each CRISPRi plate was subjected to two different condition
(Basal and 150 mM of Acetic acid). Therefore, for each plate four
different files are generated. All files generated in a single
independent experimental round are stored in a single folder.

  - **SOM\_SCR\_R001** : Raw data for round1

  - **SOM\_SCR\_R002** : Raw data for round2

**ABSOLUTE DATA**

The **Absolute** dataset gives the extracted phenotypes without any
spatial normalization

**NORMALIZED DATA**

The **Normalized** dataset is generated after removal of any spatial
bias. This is in log2 scale and referred as Log Strain Coefficient (LSC)
values

**FILE NAMING**

Each file is named with the plate identifier in such a way so that it
can be easily called programmatically

Eg. Plate 1 absolute data in basal (Ctrl) condition have the following
string  
*Ctrl1.phenotypes.Absolute*  
AND  
Plate 1 Normalized data in acetic acid (aa) stress have the string  
*aa1.phenotypes.Normalized*

## PURPOSE

At the end of this data import session, a single data.frame will be
generated with the data of 24 plates. The whole dataset will be labeled
with the strains attributes using the metadata key file (provided in the
COMPILED\_DATA folder)

**METADATA KEY FILE** : library\_keyfile1536.csv

#### IMPORTING THE METADATA FILE

``` r
Metadata_CRISPRi <- read.csv("COMPILED_DATA/library_keyfile1536.csv", na.strings = "#N/A", stringsAsFactors = FALSE)
str(Metadata_CRISPRi)
```

    ## 'data.frame':    36864 obs. of  11 variables:
    ##  $ SL_No             : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ gRNA_name         : chr  "FBP26-TRg-1" "FBP26-TRg-1" "HMI1-NRg-1" "HMI1-NRg-1" ...
    ##  $ Seq               : chr  "GCTTATCATACATTTACATC" "GCTTATCATACATTTACATC" "AAAAATTCTGACACATCACA" "AAAAATTCTGACACATCACA" ...
    ##  $ SOURCEPLATEID     : chr  "R2877.H.001" "R2877.H.001" "R2877.H.001" "R2877.H.001" ...
    ##  $ SOURCEDENSITY     : chr  "384A" "384A" "384A" "384A" ...
    ##  $ SOURCECOLONYCOLUMN: int  1 1 2 2 3 3 4 4 5 5 ...
    ##  $ SOURCECOLONYROW   : chr  "A" "A" "A" "A" ...
    ##  $ border            : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ GENE              : chr  "FBP26" "FBP26" "HMI1" "HMI1" ...
    ##  $ Control.gRNA      : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536     : chr  "A1" "A2" "A3" "A4" ...

#### GENERATE BASAL **ABSOLUTE** DATASET

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_Ctrl_Abs <- data.frame()
for(i in 1:24){
  m <- paste0("Ctrl", i, ".phenotypes.Absolute") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_Ctrl_Abs <- rbind(data_Ctrl_Abs, temp_df)
}
str(data_Ctrl_Abs)
```

    ## 'data.frame':    36864 obs. of  18 variables:
    ##  $ Plate                                   : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Row                                     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Column                                  : int  0 1 2 3 4 5 6 7 8 9 ...
    ##  $ Phenotypes.InitialValue                 : num  80116 83495 92467 97104 99622 ...
    ##  $ Phenotypes.ExperimentBaseLine           : num  81419 84504 94039 98641 101354 ...
    ##  $ Phenotypes.ExperimentEndAverage         : num  8615209 6710686 6037319 5554230 5234758 ...
    ##  $ Phenotypes.ColonySize48h                : num  6487316 5329111 4946774 4684166 4440201 ...
    ##  $ Phenotypes.ChapmanRichardsParam2        : num  14.5 12.8 42.7 17.6 18.3 ...
    ##  $ Phenotypes.ChapmanRichardsParam3        : num  -2.74 -2.65 -2.58 -2.52 -2.47 ...
    ##  $ Phenotypes.ChapmanRichardsParamXtra     : num  16 16 16.2 16.3 16.3 ...
    ##  $ Phenotypes.ChapmanRichardsParam1        : num  1.95 1.89 1.84 1.8 1.77 ...
    ##  $ Phenotypes.ChapmanRichardsParam4        : num  -31.33 -31.47 -3.93 -2.78 -2.48 ...
    ##  $ Phenotypes.GenerationTimeStErrOfEstimate: num  0.012879 0.002174 0.000633 0.001667 0.000759 ...
    ##  $ Phenotypes.ExperimentGrowthYield        : num  8533790 6626182 5943280 5455589 5133404 ...
    ##  $ Phenotypes.GenerationTime               : num  2.53 2.48 2.51 2.56 2.53 ...
    ##  $ Phenotypes.ExperimentPopulationDoublings: num  6.73 6.31 6 5.82 5.69 ...
    ##  $ Phenotypes.ChapmanRichardsFit           : num  0.999 0.999 0.999 0.999 0.999 ...
    ##  $ Phenotypes.GenerationTimeWhen           : num  6.14 4.1 4.1 3.76 3.76 ...
