CRISPRi PHENOMICS DATA ANALYSIS
================
Vaskar Mukherjee
2/3/2021

  - [IMPORT SCAN-O-MATIC RAW DATA](#import-scan-o-matic-raw-data)
      - [PURPOSE1](#purpose1)
      - [IMPORTING THE METADATA FILE](#importing-the-metadata-file)
      - [GENERATE BASAL **ABSOLUTE**
        DATASET](#generate-basal-absolute-dataset)
      - [GENERATE ACETIC ACID **ABSOLUTE**
        DATASET](#generate-acetic-acid-absolute-dataset)
      - [GENERATE BASAL **NORMALIZED**
        DATASET](#generate-basal-normalized-dataset)
      - [GENERATE ACETIC ACID **NORMALIZED**
        DATASET](#generate-acetic-acid-normalized-dataset)
      - [COMBINE THE DATASETS TO OBTAIN FINAL
        DATAFRAME](#combine-the-datasets-to-obtain-final-dataframe)
      - [IMPORT RESULTS FROM ROUND1](#import-results-from-round1)
      - [COMBINE THE DATASETS of ROUND 1 AND
        2](#combine-the-datasets-of-round-1-and-2)
  - [SCAN-O-MATIC PHENOMICS ANALYSIS](#scan-o-matic-phenomics-analysis)
      - [PURPOSE2](#purpose2)
      - [ESTIMATE THE LOG PHENOTYPIC INDEX (LPI)
        VALUES](#estimate-the-log-phenotypic-index-lpi-values)
      - [PERFORM PLATE-WISE BATCH
        CORRECTION](#perform-plate-wise-batch-correction)
      - [ESTIMATE THE BATCH CORRECTED LOG PHENOTYPIC INDEX (LPI)
        VALUES](#estimate-the-batch-corrected-log-phenotypic-index-lpi-values)
      - [SETTING THE NAMES OF THE NEW
        COLUMNS](#setting-the-names-of-the-new-columns)
      - [EXTRACT ONLY THE BATCH CORRECTED
        COLUMNS](#extract-only-the-batch-corrected-columns)
      - [CONSTRUCT A NEW DATA
        STRUCTURE](#construct-a-new-data-structure)
          - [REMOVE ROWS WITH SPATIAL CONTROL STRAIN
            DATA](#remove-rows-with-spatial-control-strain-data)
          - [CREATE A TABLE OF UNIQUE
            gRNA](#create-a-table-of-unique-grna)
          - [ARRANGE THE DATA IN THE DESIRED
            FORMAT](#arrange-the-data-in-the-desired-format)
          - [ASSIGN COLUMN NAMES](#assign-column-names)
      - [PERFORM STATISTICAL ANALYSIS](#perform-statistical-analysis)
          - [METHOD 1](#method-1)
              - [RECALCULATION OF SOME PHENOTYPIC
                PARAMETERS](#recalculation-of-some-phenotypic-parameters)
              - [EXTRACT ALL LPI GT MEAN DATA POINTS WITHIN
                INTER-QUARTILE-RANGE
                (IQR)](#extract-all-lpi-gt-mean-data-points-within-inter-quartile-range-iqr)
              - [ESTIMATE P-VALUE](#estimate-p-value)
              - [FALSE DISCOVERY RATE ADJUSTMENT OF
                P-VALUE](#false-discovery-rate-adjustment-of-p-value)
              - [P-VALUE DISGNOSTICS FOR
                METHOD1](#p-value-disgnostics-for-method1)

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

## PURPOSE1

At the end of this data import session, a single data.frame will be
generated with the data of 24 plates. The whole dataset will be labeled
with the strains attributes using the metadata key file (provided in the
COMPILED\_DATA folder). The data import below is shown for only Round2
dataset. Round1 can be generated modifying the folder location

**METADATA KEY FILE** : library\_keyfile1536.csv

## IMPORTING THE METADATA FILE

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

## GENERATE BASAL **ABSOLUTE** DATASET

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

Several phenotypes are extracted. However, the most useful for this
study will be,

  - Column No: 14 i.e. **Phenotypes.ExperimentGrowthYield**
  - Column No: 15 i.e. **Phenotypes.GenerationTime**

Extract only this two column in the final data.frame  
Rename the column names to prevent any ambiguity

``` r
data_Ctrl_Abs_Trim <- data_Ctrl_Abs[, 14:15]
colnames(data_Ctrl_Abs_Trim) <- c("CTRL_Y_ABS", "CTRL_GT_ABS")
str(data_Ctrl_Abs_Trim)
```

    ## 'data.frame':    36864 obs. of  2 variables:
    ##  $ CTRL_Y_ABS : num  8533790 6626182 5943280 5455589 5133404 ...
    ##  $ CTRL_GT_ABS: num  2.53 2.48 2.51 2.56 2.53 ...

## GENERATE ACETIC ACID **ABSOLUTE** DATASET

Following the same strategy as above

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_AA_Abs <- data.frame()
for(i in 1:24){
  m <- paste0("aa", i, ".phenotypes.Absolute") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_AA_Abs <- rbind(data_AA_Abs, temp_df)
}
data_AA_Abs_Trim <- data_AA_Abs[, 14:15]
colnames(data_AA_Abs_Trim) <- c("AA_Y_ABS", "AA_GT_ABS")
```

## GENERATE BASAL **NORMALIZED** DATASET

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_Ctrl_Norm <- data.frame()
for(i in 1:24){
  m <- paste0("Ctrl", i, ".phenotypes.Normalized") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_Ctrl_Norm <- rbind(data_Ctrl_Norm, temp_df)
}
str(data_Ctrl_Norm)
```

    ## 'data.frame':    36864 obs. of  8 variables:
    ##  $ Plate                                   : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Row                                     : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Column                                  : int  0 1 2 3 4 5 6 7 8 9 ...
    ##  $ Phenotypes.ExperimentGrowthYield        : num  0.755 0.39 0.233 0.484 0.396 ...
    ##  $ Phenotypes.GenerationTime               : num  0.0505 0.0204 0.0389 -0.0173 -0.0327 ...
    ##  $ Phenotypes.ExperimentPopulationDoublings: num  0.1907 0.0991 0.0272 0.113 0.0818 ...
    ##  $ Phenotypes.ExperimentBaseLine           : num  -0.0884 -0.0347 0.1195 0.0367 0.0758 ...
    ##  $ Phenotypes.ColonySize48h                : num  0.584 0.301 0.193 0.397 0.32 ...

The most useful for this study will be,

  - Column No: 4 i.e. Phenotypes.ExperimentGrowthYield
  - Column No: 5 i.e. Phenotypes.GenerationTime

Extract only this two column

``` r
data_Ctrl_Norm_Trim <- data_Ctrl_Norm[, 4:5]
colnames(data_Ctrl_Norm_Trim) <- c("CTRL_Y_NORM", "CTRL_GT_NORM")
```

## GENERATE ACETIC ACID **NORMALIZED** DATASET

Same as above

``` r
m <- vector(mode = "character", length = 0)
file.names<-vector(mode = "character", length = 0)
temp_df<-data.frame()
data_AA_Norm <- data.frame()
for(i in 1:24){
  m <- paste0("aa", i, ".phenotypes.Normalized") 
  file.names[i] <- dir("RAW_DATA/SOM_SCR_R002/", pattern = m, full.names = TRUE)
  temp_df <- read.csv(file.names[i], na.strings = "NoGrowth")
  data_AA_Norm <- rbind(data_AA_Norm, temp_df)
}
data_AA_Norm_Trim <- data_AA_Norm[, 4:5]
colnames(data_AA_Norm_Trim) <- c("AA_Y_NORM", "AA_GT_NORM")
```

## COMBINE THE DATASETS TO OBTAIN FINAL DATAFRAME

Trimmed datasets are combined to obtain the final data.frame. The
combined data frame is labeled as data from ROUND2

``` r
R <- rep("2nd_round", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_R2 <- cbind(Metadata_CRISPRi, 
                       Round_ID, 
                       data_Ctrl_Abs_Trim, 
                       data_AA_Abs_Trim, 
                       data_Ctrl_Norm_Trim, 
                       data_AA_Norm_Trim)
colnames(whole_data_R2)[12] <- "Round_ID"
str(whole_data_R2)
```

    ## 'data.frame':    36864 obs. of  20 variables:
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
    ##  $ Round_ID          : chr  "2nd_round" "2nd_round" "2nd_round" "2nd_round" ...
    ##  $ CTRL_Y_ABS        : num  8533790 6626182 5943280 5455589 5133404 ...
    ##  $ CTRL_GT_ABS       : num  2.53 2.48 2.51 2.56 2.53 ...
    ##  $ AA_Y_ABS          : num  2090439 2241914 1861277 1920070 1957912 ...
    ##  $ AA_GT_ABS         : num  8.92 8.43 9.19 8.91 8.87 ...
    ##  $ CTRL_Y_NORM       : num  0.755 0.39 0.233 0.484 0.396 ...
    ##  $ CTRL_GT_NORM      : num  0.0505 0.0204 0.0389 -0.0173 -0.0327 ...
    ##  $ AA_Y_NORM         : num  0.373 0.474 0.205 0.386 0.415 ...
    ##  $ AA_GT_NORM        : num  0.019 -0.0614 0.0621 -0.2135 -0.219 ...

## IMPORT RESULTS FROM ROUND1

The results from Round1 is already compiled to a .csv file in
COMPILED\_DATA folder **Results 1st Round** :
20190903\_CRISPRi\_Screen\_aa\_1st\_round.csv

Import the dataset and label as data from ROUND1

``` r
First_round <- read.csv("COMPILED_DATA/20190903_CRISPRi_Screen_aa_1st_round.csv", 
                        na.strings = c("#N/A", "NoGrowth"), 
                        stringsAsFactors = FALSE)
R <- rep("1st_round", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_R1 <- cbind(Metadata_CRISPRi, Round_ID, First_round[, 12:19])
colnames(whole_data_R1)[12] <- "Round_ID"
str(whole_data_R1)
```

    ## 'data.frame':    36864 obs. of  20 variables:
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
    ##  $ Round_ID          : chr  "1st_round" "1st_round" "1st_round" "1st_round" ...
    ##  $ CTRL_Y_ABS        : num  7046465 5380541 4841666 4517281 4209174 ...
    ##  $ CTRL_GT_ABS       : num  3.15 3.64 3.18 3.04 3.13 ...
    ##  $ AA_Y_ABS          : num  833657 844666 708042 734725 717147 ...
    ##  $ AA_GT_ABS         : num  13.1 13.6 14.5 14.7 13.6 ...
    ##  $ CTRL_Y_NORM       : num  0.899 0.51 0.358 0.565 0.463 ...
    ##  $ CTRL_GT_NORM      : num  -0.10076 0.10855 -0.08699 0.00504 0.0462 ...
    ##  $ AA_Y_NORM         : num  -0.362 -0.343 -0.598 -0.642 -0.677 ...
    ##  $ AA_GT_NORM        : num  0.273 0.331 0.423 0.425 0.313 ...

## COMBINE THE DATASETS of ROUND 1 AND 2

``` r
whole_data_CRISPRi_aa <- rbind(whole_data_R1, whole_data_R2)
```

# SCAN-O-MATIC PHENOMICS ANALYSIS

In this study most of the downstream analysis was performed using the
phenotype Generation\_time(GT)

## PURPOSE2

Downstream data processing and statistical analysis of SCAN-O-MATIC raw
output

## ESTIMATE THE LOG PHENOTYPIC INDEX (LPI) VALUES

LPI of strain is the difference of its normalized Generation\_Time(GT) /
Yield(Y) (LSC, see [IMPORT SCAN-O-MATIC RAW
DATA](#import-scan-o-matic-raw-data)) on acetic acid stress plate to the
basal condition. It gives a **RELATIVE** estimate of how a strain
performed under acetic acid stress relative to the basal condition.

The **RELATIVE GENERATION TIME** i.e. LPI\_GT = LSC\_GT\_Acetic\_Acid -
LSC\_GT\_Basal

``` r
whole_data_CRISPRi_aa[, 21] <- whole_data_CRISPRi_aa[, 19]-whole_data_CRISPRi_aa[, 17]
whole_data_CRISPRi_aa[, 22] <- whole_data_CRISPRi_aa[, 20]-whole_data_CRISPRi_aa[, 18]
colnames(whole_data_CRISPRi_aa)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_aa)[22] <- "LPI_GT"
```

## PERFORM PLATE-WISE BATCH CORRECTION

Plate-wise batch correction was conducted by subtracting the median of
LSC GT values of all the individual colonies on a plate from the
individual LSC GT values of the colonies growing on that plate.

i.e. if strainX is growing in Basal condition on plate Z, the corrected
LSC\_GT value for strainX in the Basal condition is the following;

  - LSC\_GT\_Basal\_Corrected<sub>strainX</sub> =
    (LSC\_GT\_Basal<sub>strainX</sub>) - Median(LSC\_GT
    Basal<sub>PlateZ</sub>)

<!-- end list -->

``` r
plate_ID <- as.character(unique(whole_data_CRISPRi_aa$SOURCEPLATEID))
whole_data_CRISPRi_aa_corrected <- whole_data_CRISPRi_aa
med_LogLSCctrl_RND1_GT <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND1_GT <- vector(mode = "integer", length = 0)
med_LogLSCctrl_RND2_GT <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND2_GT <- vector(mode = "integer", length = 0)
med_LogLSCctrl_RND1_Y <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND1_Y <- vector(mode = "integer", length = 0)
med_LogLSCctrl_RND2_Y <- vector(mode = "integer", length = 0)
med_LogLSCaa_RND2_Y <- vector(mode = "integer", length = 0)

for(i in 1:24){
med_LogLSCctrl_RND1_GT[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                       & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM) 
                                                                                         & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCaa_RND1_GT[i] <- median(whole_data_CRISPRi_aa_corrected$AA_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                     & !is.na(whole_data_CRISPRi_aa_corrected$AA_GT_NORM) 
                                                                                     & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCctrl_RND2_GT[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                         & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_GT_NORM) 
                                                                                         & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])

med_LogLSCaa_RND2_GT[i] <- median(whole_data_CRISPRi_aa_corrected$AA_GT_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                     & !is.na(whole_data_CRISPRi_aa_corrected$AA_GT_NORM) 
                                                                                     & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])

med_LogLSCctrl_RND1_Y[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                       & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM) 
                                                                                       & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCaa_RND1_Y[i] <- median(whole_data_CRISPRi_aa_corrected$AA_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                   & !is.na(whole_data_CRISPRi_aa_corrected$AA_Y_NORM) 
                                                                                   & whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round")])

med_LogLSCctrl_RND2_Y[i] <- median(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                       & !is.na(whole_data_CRISPRi_aa_corrected$CTRL_Y_NORM) 
                                                                                       & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])

med_LogLSCaa_RND2_Y[i] <- median(whole_data_CRISPRi_aa_corrected$AA_Y_NORM[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i]
                                                                                   & !is.na(whole_data_CRISPRi_aa_corrected$AA_Y_NORM) 
                                                                                   & whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round")])
  
whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 23] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 17] - med_LogLSCctrl_RND1_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 24] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 18] - med_LogLSCctrl_RND1_GT[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 23] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 17] - med_LogLSCctrl_RND2_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 24] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 18] - med_LogLSCctrl_RND2_GT[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 25] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 19] - med_LogLSCaa_RND1_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round") , 26] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="1st_round"), 20] - med_LogLSCaa_RND1_GT[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 25] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 19] - med_LogLSCaa_RND2_Y[i]
  whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                          whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round") , 26] <- whole_data_CRISPRi_aa_corrected[which(whole_data_CRISPRi_aa_corrected$SOURCEPLATEID==plate_ID[i] &
                                                                                                                                                  whole_data_CRISPRi_aa_corrected$Round_ID=="2nd_round"), 20] - med_LogLSCaa_RND2_GT[i]
}
```

## ESTIMATE THE BATCH CORRECTED LOG PHENOTYPIC INDEX (LPI) VALUES

Estimate the corrected LPI values (see [ESTIMATE THE LOG PHENOTYPIC
INDEX (LPI) VALUES](#estimate-the-log-phenotypic-index-lpi-values))
based on the corrected LSC values

i.e. LPI\_GT<sub>corrected</sub> =
LSC\_GT\_Acetic\_Acid<sub>corrected</sub> -
LSC\_GT\_Basal<sub>corrected</sub>

**Estimate the corrected LPI\_Y**

``` r
whole_data_CRISPRi_aa_corrected[, 27] <- whole_data_CRISPRi_aa_corrected[, 25] - whole_data_CRISPRi_aa_corrected[, 23]
```

**Estimate the corrected LPI\_GT**

``` r
whole_data_CRISPRi_aa_corrected[, 28] <- whole_data_CRISPRi_aa_corrected[, 26] - whole_data_CRISPRi_aa_corrected[, 24] 
```

## SETTING THE NAMES OF THE NEW COLUMNS

``` r
colnm <- colnames(whole_data_CRISPRi_aa)[17:22]
colnm <- paste0(colnm, "_CR")
colnames(whole_data_CRISPRi_aa_corrected)[23:28] <- colnm
str(whole_data_CRISPRi_aa_corrected)
```

    ## 'data.frame':    73728 obs. of  28 variables:
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
    ##  $ Round_ID          : chr  "1st_round" "1st_round" "1st_round" "1st_round" ...
    ##  $ CTRL_Y_ABS        : num  7046465 5380541 4841666 4517281 4209174 ...
    ##  $ CTRL_GT_ABS       : num  3.15 3.64 3.18 3.04 3.13 ...
    ##  $ AA_Y_ABS          : num  833657 844666 708042 734725 717147 ...
    ##  $ AA_GT_ABS         : num  13.1 13.6 14.5 14.7 13.6 ...
    ##  $ CTRL_Y_NORM       : num  0.899 0.51 0.358 0.565 0.463 ...
    ##  $ CTRL_GT_NORM      : num  -0.10076 0.10855 -0.08699 0.00504 0.0462 ...
    ##  $ AA_Y_NORM         : num  -0.362 -0.343 -0.598 -0.642 -0.677 ...
    ##  $ AA_GT_NORM        : num  0.273 0.331 0.423 0.425 0.313 ...
    ##  $ LPI_Y             : num  -1.262 -0.854 -0.956 -1.207 -1.14 ...
    ##  $ LPI_GT            : num  0.374 0.223 0.51 0.42 0.266 ...
    ##  $ CTRL_Y_NORM_CR    : num  0.918 0.529 0.376 0.583 0.482 ...
    ##  $ CTRL_GT_NORM_CR   : num  -0.0877 0.1216 -0.074 0.0181 0.0592 ...
    ##  $ AA_Y_NORM_CR      : num  -0.0736 -0.0547 -0.3092 -0.3534 -0.3883 ...
    ##  $ AA_GT_NORM_CR     : num  0.19 0.248 0.34 0.342 0.229 ...
    ##  $ LPI_Y_CR          : num  -0.991 -0.583 -0.686 -0.937 -0.87 ...
    ##  $ LPI_GT_CR         : num  0.278 0.127 0.414 0.323 0.17 ...

## EXTRACT ONLY THE BATCH CORRECTED COLUMNS

``` r
whole_data_CRISPRi_aa_2 <- whole_data_CRISPRi_aa_corrected[, c(1:16, 23:28)]
colnames(whole_data_CRISPRi_aa_2)[17:22] <- colnames(whole_data_CRISPRi_aa)[17:22]
str(whole_data_CRISPRi_aa_2)
```

    ## 'data.frame':    73728 obs. of  22 variables:
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
    ##  $ Round_ID          : chr  "1st_round" "1st_round" "1st_round" "1st_round" ...
    ##  $ CTRL_Y_ABS        : num  7046465 5380541 4841666 4517281 4209174 ...
    ##  $ CTRL_GT_ABS       : num  3.15 3.64 3.18 3.04 3.13 ...
    ##  $ AA_Y_ABS          : num  833657 844666 708042 734725 717147 ...
    ##  $ AA_GT_ABS         : num  13.1 13.6 14.5 14.7 13.6 ...
    ##  $ CTRL_Y_NORM       : num  0.918 0.529 0.376 0.583 0.482 ...
    ##  $ CTRL_GT_NORM      : num  -0.0877 0.1216 -0.074 0.0181 0.0592 ...
    ##  $ AA_Y_NORM         : num  -0.0736 -0.0547 -0.3092 -0.3534 -0.3883 ...
    ##  $ AA_GT_NORM        : num  0.19 0.248 0.34 0.342 0.229 ...
    ##  $ LPI_Y             : num  -0.991 -0.583 -0.686 -0.937 -0.87 ...
    ##  $ LPI_GT            : num  0.278 0.127 0.414 0.323 0.17 ...

## CONSTRUCT A NEW DATA STRUCTURE

Construct a new data structure where data from each strain (have a
unique guide-RNA) is in a separate row and the replicates from first and
second round are side by side. Add also the mean, median and standard
deviation statistics for each phenotype

### REMOVE ROWS WITH SPATIAL CONTROL STRAIN DATA

``` r
Data_CRISPRi_aa <- subset(whole_data_CRISPRi_aa_2, whole_data_CRISPRi_aa_2$gRNA_name!="SP_Ctrl_CC23")
```

### CREATE A TABLE OF UNIQUE gRNA

``` r
df_unique_sgRNA <- data.frame(table(Data_CRISPRi_aa$gRNA_name))
```

### ARRANGE THE DATA IN THE DESIRED FORMAT

``` r
R1<-vector(mode = "integer", length = 0)
R2<-vector(mode = "integer", length = 0)
test2<-data.frame()
n<-nrow(df_unique_sgRNA)
for(i in 1:n){
  R1 <- which(Data_CRISPRi_aa$gRNA_name==df_unique_sgRNA$Var1[i] & Data_CRISPRi_aa$Round_ID=="1st_round")
  R2 <- which(Data_CRISPRi_aa$gRNA_name==df_unique_sgRNA$Var1[i] & Data_CRISPRi_aa$Round_ID=="2nd_round")
  test1 <- Data_CRISPRi_aa[c(R1, R2), ]
  test2[i, c(1:8)]<-test1[1, c(2:4, 6:7, 9:11)]
  test2[i, c(9:14)] <- test1$CTRL_GT_NORM
  test2[i, 15] <- mean(test1$CTRL_GT_NORM[1:3])
  test2[i, 16] <- mean(test1$CTRL_GT_NORM[4:6])
  test2[i, 17] <- sd(test1$CTRL_GT_NORM[1:3])
  test2[i, 18] <- sd(test1$CTRL_GT_NORM[4:6])
  test2[i, 19] <- mean(test1$CTRL_GT_NORM[1:6])
  test2[i, 20] <- median(test1$CTRL_GT_NORM[1:6])
  test2[i, 21] <- sd(test1$CTRL_GT_NORM[1:6])
  test2[i, c(22:27)] <- test1$AA_GT_NORM
  test2[i, 28] <- mean(test1$AA_GT_NORM[1:3])
  test2[i, 29] <- mean(test1$AA_GT_NORM[4:6])
  test2[i, 30] <- sd(test1$AA_GT_NORM[1:3])
  test2[i, 31] <- sd(test1$AA_GT_NORM[4:6])
  test2[i, 32] <- mean(test1$AA_GT_NORM[1:6])
  test2[i, 33] <- median(test1$AA_GT_NORM[1:6])
  test2[i, 34] <- sd(test1$AA_GT_NORM[1:6])
  test2[i, c(35:40)] <- test1$LPI_GT
  test2[i, 41] <- mean(test1$LPI_GT[1:3])
  test2[i, 42] <- mean(test1$LPI_GT[4:6])
  test2[i, 43] <- sd(test1$LPI_GT[1:3])
  test2[i, 44] <- sd(test1$LPI_GT[4:6])
  test2[i, 45] <- mean(test1$LPI_GT[1:6])
  test2[i, 46] <- median(test1$LPI_GT[1:6])
  test2[i, 47] <- sd(test1$LPI_GT[1:6])
  test2[i, c(48:53)] <- test1$CTRL_Y_NORM
  test2[i, 54] <- mean(test1$CTRL_Y_NORM[1:3])
  test2[i, 55] <- mean(test1$CTRL_Y_NORM[4:6])
  test2[i, 56] <- sd(test1$CTRL_Y_NORM[1:3])
  test2[i, 57] <- sd(test1$CTRL_Y_NORM[4:6])
  test2[i, 58] <- mean(test1$CTRL_Y_NORM[1:6])
  test2[i, 59] <- median(test1$CTRL_Y_NORM[1:6])
  test2[i, 60] <- sd(test1$CTRL_Y_NORM[1:6])
  test2[i, c(61:66)] <- test1$AA_Y_NORM
  test2[i, 67] <- mean(test1$AA_Y_NORM[1:3])
  test2[i, 68] <- mean(test1$AA_Y_NORM[4:6])
  test2[i, 69] <- sd(test1$AA_Y_NORM[1:3])
  test2[i, 70] <- sd(test1$AA_Y_NORM[4:6])
  test2[i, 71] <- mean(test1$AA_Y_NORM[1:6])
  test2[i, 72] <- median(test1$AA_Y_NORM[1:6])
  test2[i, 73] <- sd(test1$AA_Y_NORM[1:6])
  test2[i, c(74:79)] <- test1$LPI_Y
  test2[i, 80] <- mean(test1$LPI_Y[1:3])
  test2[i, 81] <- mean(test1$LPI_Y[4:6])
  test2[i, 82] <- sd(test1$LPI_Y[1:3])
  test2[i, 83] <- sd(test1$LPI_Y[4:6])
  test2[i, 84] <- mean(test1$LPI_Y[1:6])
  test2[i, 85] <- median(test1$LPI_Y[1:6])
  test2[i, 86] <- sd(test1$LPI_Y[1:6])
}
```

### ASSIGN COLUMN NAMES

Column names are already stored in a text times available in the
**COMPILED\_DATA** folder. Then store the data.frame under a new name.

``` r
column_names <- read.table("COMPILED_DATA/Column_names.txt", header = FALSE, sep = "\t", as.is = TRUE)
colnames(test2) <- column_names$V1
Analysis_CRISPRi_aa_Complete <- test2
str(Analysis_CRISPRi_aa_Complete)
```

    ## 'data.frame':    9078 obs. of  86 variables:
    ##  $ gRNA_name          : chr  "AAR2-NRg-3" "AAR2-NRg-4" "AAR2-TRg-15" "AAR2-TRg-16" ...
    ##  $ Seq                : chr  "CCAGCGATAAGGAGGATCTT" "TGTGTCCTTTCTTCATCTCT" "AAAAGGAAAAAGTAATTAGG" "GTGAAAAGGAAAAAGTAATT" ...
    ##  $ SOURCEPLATEID      : chr  "R2877.H.023" "R2877.H.024" "R2877.H.023" "R2877.H.023" ...
    ##  $ SOURCECOLONYCOLUMN : int  5 21 22 20 21 18 6 7 9 8 ...
    ##  $ SOURCECOLONYROW    : chr  "O" "L" "P" "N" ...
    ##  $ GENE               : chr  "AAR2" "AAR2" "AAR2" "AAR2" ...
    ##  $ Control.gRNA       : int  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ Location_1536      : chr  "AC9" "W41" "AE43" "AA39" ...
    ##  $ CTRL_GT_RND1_R1    : num  0.00397 0.00586 0.10456 -0.01142 0.0174 ...
    ##  $ CTRL_GT_RND1_R2    : num  -0.0086 0.00166 0.07118 0.00256 0.04094 ...
    ##  $ CTRL_GT_RND1_R3    : num  -0.000625 -0.043918 -0.009614 0.007969 -0.014018 ...
    ##  $ CTRL_GT_RND2_R1    : num  0.0431 -0.0336 0.0292 0.0381 0.0742 ...
    ##  $ CTRL_GT_RND2_R2    : num  -0.0165 0.00744 0.01433 0.06599 -0.02973 ...
    ##  $ CTRL_GT_RND2_R3    : num  0.03266 -0.04468 -0.04455 0.01994 -0.00924 ...
    ##  $ CTRL_GT_RND1_MEAN  : num  -0.00175 -0.01213 0.05538 -0.0003 0.01477 ...
    ##  $ CTRL_GT_RND2_MEAN  : num  0.019756 -0.023599 -0.000341 0.041342 0.011747 ...
    ##  $ CTRL_GT_RND1_SD    : num  0.00636 0.02761 0.05871 0.01001 0.02757 ...
    ##  $ CTRL_GT_RND2_SD    : num  0.0318 0.0275 0.039 0.0232 0.0551 ...
    ##  $ CTRL_GT_RND1_2_MEAN: num  0.009 -0.0179 0.0275 0.0205 0.0133 ...
    ##  $ CTRL_GT_RND1_2_MED : num  0.00167 -0.01595 0.02176 0.01395 0.00408 ...
    ##  $ CTRL_GT_RND1_2_SD  : num  0.0237 0.0254 0.054 0.0278 0.039 ...
    ##  $ AA_GT_RND1_R1      : num  -0.0239 0.06 -0.0178 -0.1839 0.5993 ...
    ##  $ AA_GT_RND1_R2      : num  -0.0311 0.0763 -0.0145 -0.1293 0.2619 ...
    ##  $ AA_GT_RND1_R3      : num  0.0315 -0.0191 -0.0778 -0.0514 0.3276 ...
    ##  $ AA_GT_RND2_R1      : num  0.0142 -0.1375 0.111 0.0105 0.2021 ...
    ##  $ AA_GT_RND2_R2      : num  0.0383 -0.2266 0.0588 0.0966 0.1831 ...
    ##  $ AA_GT_RND2_R3      : num  0.0265 -0.2298 0.0789 0.0787 0.0114 ...
    ##  $ AA_GT_RND1_MEAN    : num  -0.00781 0.03905 -0.03671 -0.12151 0.39624 ...
    ##  $ AA_GT_RND2_MEAN    : num  0.0263 -0.198 0.0829 0.0619 0.1322 ...
    ##  $ AA_GT_RND1_SD      : num  0.0343 0.051 0.0356 0.0666 0.1789 ...
    ##  $ AA_GT_RND2_SD      : num  0.0121 0.0524 0.0263 0.0455 0.105 ...
    ##  $ AA_GT_RND1_2_MEAN  : num  0.00925 -0.07945 0.02308 -0.02979 0.2642 ...
    ##  $ AA_GT_RND1_2_MED   : num  0.0203 -0.0783 0.0221 -0.0205 0.232 ...
    ##  $ AA_GT_RND1_2_SD    : num  0.0296 0.1378 0.0712 0.1127 0.1953 ...
    ##  $ LPI_GT_RND1_R1     : num  -0.0278 0.0541 -0.1224 -0.1724 0.5819 ...
    ##  $ LPI_GT_RND1_R2     : num  -0.0225 0.0746 -0.0857 -0.1318 0.2209 ...
    ##  $ LPI_GT_RND1_R3     : num  0.0322 0.0248 -0.0682 -0.0594 0.3416 ...
    ##  $ LPI_GT_RND2_R1     : num  -0.0289 -0.1039 0.0818 -0.0276 0.1279 ...
    ##  $ LPI_GT_RND2_R2     : num  0.0548 -0.2341 0.0444 0.0306 0.2128 ...
    ##  $ LPI_GT_RND2_R3     : num  -0.0062 -0.1851 0.1234 0.0588 0.0206 ...
    ##  $ LPI_GT_RND1_MEAN   : num  -0.00605 0.05119 -0.09208 -0.12121 0.38147 ...
    ##  $ LPI_GT_RND2_MEAN   : num  0.00656 -0.17435 0.08321 0.02058 0.12042 ...
    ##  $ LPI_GT_RND1_SD     : num  0.0332 0.025 0.0276 0.0573 0.1837 ...
    ##  $ LPI_GT_RND2_SD     : num  0.0433 0.0657 0.0395 0.0441 0.0963 ...
    ##  $ LPI_GT_RND1_2_MEAN : num  0.000251 -0.061582 -0.004437 -0.050313 0.250944 ...
    ##  $ LPI_GT_RND1_2_MED  : num  -0.0143 -0.0396 -0.0119 -0.0435 0.2169 ...
    ##  $ LPI_GT_RND1_2_SD   : num  0.0352 0.1313 0.1007 0.0901 0.1941 ...
    ##  $ CTRL_Y_RND1_R1     : num  0.055 0.0573 0.0189 0.0779 -0.1151 ...
    ##  $ CTRL_Y_RND1_R2     : num  0.0131 0.0472 0.0129 0.0562 -0.0723 ...
    ##  $ CTRL_Y_RND1_R3     : num  0.0306 0.0109 -0.0236 0.1208 -0.0171 ...
    ##  $ CTRL_Y_RND2_R1     : num  0.0399 0.0113 -0.1422 0.0605 -0.0676 ...
    ##  $ CTRL_Y_RND2_R2     : num  0.02531 0.0066 -0.17285 -0.00976 -0.08442 ...
    ##  $ CTRL_Y_RND2_R3     : num  -0.0528 0.0089 0.00598 0.06854 -0.08083 ...
    ##  $ CTRL_Y_RND1_MEAN   : num  0.03288 0.03847 0.00276 0.08497 -0.06817 ...
    ##  $ CTRL_Y_RND2_MEAN   : num  0.00414 0.00893 -0.10301 0.03977 -0.07762 ...
    ##  $ CTRL_Y_RND1_SD     : num  0.0211 0.0244 0.023 0.0329 0.0491 ...
    ##  $ CTRL_Y_RND2_SD     : num  0.04985 0.00234 0.09563 0.04308 0.00885 ...
    ##  $ CTRL_Y_RND1_2_MEAN : num  0.0185 0.0237 -0.0501 0.0624 -0.0729 ...
    ##  $ CTRL_Y_RND1_2_MED  : num  0.028 0.0111 -0.0088 0.0645 -0.0766 ...
    ##  $ CTRL_Y_RND1_2_SD   : num  0.0377 0.0224 0.085 0.0423 0.032 ...
    ##  $ AA_Y_RND1_R1       : num  0.0672 0.0106 0.1785 0.272 -2.1 ...
    ##  $ AA_Y_RND1_R2       : num  -0.3832 0.0102 0.1196 0.232 -1.1873 ...
    ##  $ AA_Y_RND1_R3       : num  -0.1599 -0.0465 -0.0381 0.1036 -1.0143 ...
    ##  $ AA_Y_RND2_R1       : num  0.0503 0.3083 -0.2434 0.0968 -0.4223 ...
    ##  $ AA_Y_RND2_R2       : num  -0.0505 0.4391 -0.2526 -0.094 -0.304 ...
    ##  $ AA_Y_RND2_R3       : num  -0.00342 0.52795 -0.20049 -0.02924 -0.26528 ...
    ##  $ AA_Y_RND1_MEAN     : num  -0.15864 -0.00857 0.08665 0.20251 -1.43385 ...
    ##  $ AA_Y_RND2_MEAN     : num  -0.00121 0.42511 -0.23216 -0.00879 -0.33052 ...
    ##  $ AA_Y_RND1_SD       : num  0.2252 0.0329 0.112 0.088 0.5834 ...
    ##  $ AA_Y_RND2_SD       : num  0.0504 0.1105 0.0278 0.097 0.0818 ...
    ##  $ AA_Y_RND1_2_MEAN   : num  -0.0799 0.2083 -0.0728 0.0969 -0.8822 ...
    ##  $ AA_Y_RND1_2_MED    : num  -0.027 0.159 -0.119 0.1 -0.718 ...
    ##  $ AA_Y_RND1_2_SD     : num  0.17 0.248 0.189 0.142 0.71 ...
    ##  $ LPI_Y_RND1_R1      : num  0.0122 -0.0467 0.1595 0.1942 -1.9849 ...
    ##  $ LPI_Y_RND1_R2      : num  -0.3963 -0.0369 0.1066 0.1757 -1.1149 ...
    ##  $ LPI_Y_RND1_R3      : num  -0.1905 -0.0575 -0.0145 -0.0172 -0.9972 ...
    ##  $ LPI_Y_RND2_R1      : num  0.0104 0.297 -0.1012 0.0363 -0.3547 ...
    ##  $ LPI_Y_RND2_R2      : num  -0.0758 0.4325 -0.0797 -0.0842 -0.2196 ...
    ##  $ LPI_Y_RND2_R3      : num  0.0494 0.519 -0.2065 -0.0978 -0.1844 ...
    ##  $ LPI_Y_RND1_MEAN    : num  -0.1915 -0.047 0.0839 0.1175 -1.3657 ...
    ##  $ LPI_Y_RND2_MEAN    : num  -0.00536 0.41618 -0.12915 -0.04856 -0.2529 ...
    ##  $ LPI_Y_RND1_SD      : num  0.2042 0.0103 0.0892 0.1171 0.5395 ...
    ##  $ LPI_Y_RND2_SD      : num  0.0641 0.1119 0.0678 0.0738 0.0899 ...
    ##  $ LPI_Y_RND1_2_MEAN  : num  -0.0984 0.1846 -0.0226 0.0345 -0.8093 ...
    ##  $ LPI_Y_RND1_2_MED   : num  -0.03273 0.13006 -0.04712 0.00953 -0.67593 ...
    ##  $ LPI_Y_RND1_2_SD    : num  0.169 0.263 0.137 0.126 0.701 ...

## PERFORM STATISTICAL ANALYSIS

Multiple statistical method was applied to identify the best fit
statistical model for this dataset. We start with the complete dataset
and give it a new name to avoid distorting the original dataset.

``` r
Analysis_Final <- Analysis_CRISPRi_aa_Complete
```

### METHOD 1

For **METHOD 1**, We hypothesized that the difference between the mean
phenotypic performance of a specific CRISPRi strain (StrainX) in the two
independent experimental rounds (n=2) to the mean phenotypic performance
of all the CRISPRi strains that falls within the interquartile range
(IQR) of the complete dataset would be zero, and any difference within
the IQR to be just by chance.

**Null Hypothesis** : µ(µ<sub>LPI\_GT\_StrainX\_Round1</sub>,
µ<sub>LPI\_GT\_StrainX\_Round2</sub>)- µ(InterquartileRange\_LPI\_GT) =
0

#### RECALCULATION OF SOME PHENOTYPIC PARAMETERS

In this method, we estimate the Mean / standard deviation (SD) of the
LPI GT of Round 1 and Round 2 separately for each strain. When one/two
of the three replicates of a strain in a round returned missing value
(i.e. NA), then the mean / SD of LPI GT for that round is calculated by
taking average of the non NA replicates. Therefore, excluding the
missing values the mean and SD statistics were recalculated. We
implemented a if else decision tree for this

  - The mean and SD of Normalized generation time (**LSC GT mean**) at
    **Basal condition** re-calculation

<!-- end list -->

``` r
for(i in 1:nrow(Analysis_Final)){
  x1 <- as.numeric(Analysis_Final[i, 9:11][which(!is.na(Analysis_Final[i, 9:11]))])
  x2 <- as.numeric(Analysis_Final[i, 12:14][which(!is.na(Analysis_Final[i, 12:14]))])
  if(length(x1)==0){
    Analysis_Final$CTRL_GT_RND1_MEAN[i] <- NA
  } else{
    Analysis_Final$CTRL_GT_RND1_MEAN[i] <- as.numeric(mean(x1))
  }
  if(length(x2)==0){
    Analysis_Final$CTRL_GT_RND2_MEAN[i] <- NA
  } else{
    Analysis_Final$CTRL_GT_RND2_MEAN[i] <- as.numeric(mean(x2))
  }
  if(sum(is.na(c(Analysis_Final$CTRL_GT_RND1_MEAN[i], Analysis_Final$CTRL_GT_RND2_MEAN[i])))==0){
    Analysis_Final$CTRL_GT_RND1_2_MEAN[i] <- as.numeric(mean(c(Analysis_Final$CTRL_GT_RND1_MEAN[i], Analysis_Final$CTRL_GT_RND2_MEAN[i])))
    Analysis_Final[i, 87] <- as.numeric(sd(c(Analysis_Final$CTRL_GT_RND1_MEAN[i], Analysis_Final$CTRL_GT_RND2_MEAN[i])))
  } else{
    Analysis_Final$CTRL_GT_RND1_2_MEAN[i] <- NA
    Analysis_Final[i, 87] <- NA
  }
}
colnames(Analysis_Final)[87] <- "CTRL_GT_MEAN_RND1_2_SD"
```

  - The mean and SD of Normalized generation time (**LSC GT mean**) at
    **150mM acetic acid** re-calculation

<!-- end list -->

``` r
for(i in 1:nrow(Analysis_Final)){
  x1 <- as.numeric(Analysis_Final[i, 22:24][which(!is.na(Analysis_Final[i, 22:24]))])
  x2 <- as.numeric(Analysis_Final[i, 25:27][which(!is.na(Analysis_Final[i, 25:27]))])
  if(length(x1)==0){
    Analysis_Final$AA_GT_RND1_MEAN[i] <- NA
  } else{
    Analysis_Final$AA_GT_RND1_MEAN[i] <- as.numeric(mean(x1))
  }
  if(length(x2)==0){
    Analysis_Final$AA_GT_RND2_MEAN[i] <- NA
  } else{
    Analysis_Final$AA_GT_RND2_MEAN[i] <- as.numeric(mean(x2))
  }
  if(sum(is.na(c(Analysis_Final$AA_GT_RND1_MEAN[i], Analysis_Final$AA_GT_RND2_MEAN[i])))==0){
    Analysis_Final$AA_GT_RND1_2_MEAN[i] <- as.numeric(mean(c(Analysis_Final$AA_GT_RND1_MEAN[i], Analysis_Final$AA_GT_RND2_MEAN[i])))
    Analysis_Final[i, 88] <- as.numeric(sd(c(Analysis_Final$AA_GT_RND1_MEAN[i], Analysis_Final$AA_GT_RND2_MEAN[i])))
  } else{
    Analysis_Final$AA_GT_RND1_2_MEAN[i] <- NA
    Analysis_Final[i, 88] <- NA
  }
}
colnames(Analysis_Final)[88] <- "AA_GT_MEAN_RND1_2_SD"
```

  - The mean and SD of **RELATIVE** generation time (**LPI GT mean**) at
    **150mM acetic acid** re-calculation

<!-- end list -->

``` r
for(i in 1:nrow(Analysis_Final)){
  x1 <- as.numeric(Analysis_Final[i, 35:37][which(!is.na(Analysis_Final[i, 35:37]))])
  x2 <- as.numeric(Analysis_Final[i, 38:40][which(!is.na(Analysis_Final[i, 38:40]))])
  if(length(x1)==0){
    Analysis_Final$LPI_GT_RND1_MEAN[i] <- NA
  } else{
    Analysis_Final$LPI_GT_RND1_MEAN[i] <- as.numeric(mean(x1))
  }
  if(length(x2)==0){
    Analysis_Final$LPI_GT_RND2_MEAN[i] <- NA
  } else{
    Analysis_Final$LPI_GT_RND2_MEAN[i] <- as.numeric(mean(x2))
  }
  if(sum(is.na(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))==0){
    Analysis_Final$LPI_GT_RND1_2_MEAN[i] <- as.numeric(mean(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))
    Analysis_Final[i, 89] <- as.numeric(sd(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))
  } else{
    Analysis_Final$LPI_GT_RND1_2_MEAN[i] <- NA
    Analysis_Final[i, 89] <- NA
  }
}
colnames(Analysis_Final)[89] <- "LPI_GT_MEAN_RND1_2_SD"
```

#### EXTRACT ALL LPI GT MEAN DATA POINTS WITHIN INTER-QUARTILE-RANGE (IQR)

BOX PLOT - MEAN RELATIVE GENERATION TIME (LPI GT)

![BOX PLOT - MEAN RELATIVE GENERATION TIME (LPI
GT)](CRISPRi_Screening_AceticAcid_files/figure-gfm/unnamed-chunk-25-1.png)

Display Box-plot statistics

``` r
box_stat_LPI_GT_R1_2_mean$stats
```

    ##             [,1]
    ## [1,] -0.16933911
    ## [2,] -0.02428792
    ## [3,]  0.02084505
    ## [4,]  0.07255828
    ## [5,]  0.21771505

  - 25th Percentile = -0.02428792
  - 75th Percentile = 0.07255828

Therefore, extraction of the data points within IQR

``` r
Intermediate_50 <- Analysis_Final$LPI_GT_RND1_2_MEAN[which(Analysis_Final$LPI_GT_RND1_2_MEAN>=-0.02428792
                                                           &Analysis_Final$LPI_GT_RND1_2_MEAN<=0.07255828)]
```

#### ESTIMATE P-VALUE

P-value is estimated by Welch two sample two-sided t-test (an adaptation
of Student’s t-test)

``` r
for(i in 1:nrow(Analysis_Final)){
  if(sum(is.na(c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i])))==0){
    P_value <- t.test(Intermediate_50, c(Analysis_Final$LPI_GT_RND1_MEAN[i], Analysis_Final$LPI_GT_RND2_MEAN[i]))
    Analysis_Final[i, 90] <- P_value$p.value
  } else{
    Analysis_Final[i, 90] <- NA
  }
}
colnames(Analysis_Final)[90] <- "P_value_M1"
```

#### FALSE DISCOVERY RATE ADJUSTMENT OF P-VALUE

P-value adjustment by **BENJAMINI-HOCHBERG False Discovery Rate (FDR)
method**

``` r
Analysis_Final[which(!is.na(Analysis_Final$P_value_M1)), 91] <- p.adjust(Analysis_Final$P_value_M1[which(!is.na(Analysis_Final$P_value_M1))], 
                                                                      method = "BH", 
                                                                      n = length(Analysis_Final$P_value_M1[which(!is.na(Analysis_Final$P_value_M1))]))
colnames(Analysis_Final)[91] <- "P.adjusted_M1"
```

#### P-VALUE DISGNOSTICS FOR METHOD1

NUMBER OF SIGNIFICANT STRAINS

``` r
length(Analysis_Final$P_value_M1[which(Analysis_Final$P_value_M1<=0.05)])
```

    ## [1] 434

``` r
length(Analysis_Final$P.adjusted_M1[which(Analysis_Final$P.adjusted_M1<=0.05)])
```

    ## [1] 66

``` r
length(Analysis_Final$P_value_M1[which(Analysis_Final$P_value_M1<=0.1)])
```

    ## [1] 842

``` r
length(Analysis_Final$P.adjusted_M1[which(Analysis_Final$P.adjusted_M1<=0.1)])
```

    ## [1] 71

P-VALUE DIAGNOSTICS BY **HISTOGRAM ANALYSIS**

![P-VALUE DIAGNOSTICS
METHOD1](CRISPRi_Screening_AceticAcid_files/figure-gfm/unnamed-chunk-31-1.png)
